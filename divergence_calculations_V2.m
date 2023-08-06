%The purpose of this code (after CDO process)
%1. Calculate zonal volumetric flow within the control volume and at  its borders
%2. Save zonal volumetric flows in same number of sets as input data (5)
%3. Calculate total volumetric divergence at lower limb Sigma0 boundary 
%4. Plot annually-averaged divergence and boundary flow over time 
%5. Save  divergence data and map seasonal conditions at each significant
%boundary(Labrador Sea, Boundaries around Iceland, Faroe Islands, Scotland)
clear all
close all

addpath GSW %only works if in same higher folder as GSW
addpath GSW/html
addpath GSW/library
addpath GSW/pdf
addpath GSW/thermodynamics_from_t
addpath othercolor
load colorData.mat

%% adjustable parameters
modelName = 'CESM2';
spacing = 0.5;
densitySpacing = 0.2;

%% loading .mat files
load(['regridded_dims_' modelName '_' num2str(spacing) '_V2.mat']) %loads 0.25 grid from 49N to 67N (lat), -5 to -65 (long)
load(['sigma_and_MOC_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '_V2.mat'])
%loads previous variables delta_z, depth, lat, long, time
time = []; %resetting time,not sure why

meanSigmaMOCvalue = round(mean(sigmaMocValue),2); %from sigmaMOC2023_highest_res_corrected.mat (0.25 spacing)
% meanSigmaMOCvalue = 27.66; 
[rGrid,cGrid] = size(longitude);
dy(1:rGrid,1:cGrid) = gsw_distance([1,1], [0, 0.25]); %y distance is always the same, but dx is latitude-dependent
dA(1:rGrid,1:cGrid,1:length(depth)) = 0;
D = length(depth);

for z = 1:D
    dA(:,:,z) = dy.*delta_z(z);
end

u_template = ['uo_' modelName '_regridded_subset_00000'];

sigma0_template_p1 = ['Sigma0_gridded_' modelName '_instance']; %then have instance number
sigma0_template_p2 = ['_' num2str(spacing) '_V2.nc'];
Q_template_p1 = ['Q_V_gridded_' modelName '_instance']; %then have instance number
Q_template_p2 = ['_' num2str(spacing) '_V2.nc'];

U_gridded(1:rGrid,1:cGrid,1:D) = 0;
%adding lat and long mask for small Gulf of St. Lawrence Correction
latMask = latitude >= 49 & latitude <= 52.25;
longMask = longitude>= -65 & longitude <= -56;
combinedMask = latMask & longMask;
combinedMask = repmat(combinedMask, [1, 1, D]);
multiplicationMatrix(1:rGrid,1:cGrid,1:length(depth)) = 1;
multiplicationMatrix(combinedMask)=NaN;

latMask = [];
longMask = [];
Flow_50N = [];
Flow_LS_WB = [];
Flow_LS_NB = [];
Flow_GI_NB = [];
Flow_IF_NB = [];
Flow_IF_EB = [];
Flow_GS_DB = [];

sigma0_compare_GI_NB = [];
Q_V_compare_GI_NB = [];
sigma0_compare_LS_NB = [];
Q_V_compare_LS_NB = [];
sigma0_compare_IF_NB = [];
Q_V_compare_IF_NB = [];
sigma0_compare_IF_EB = [];
Q_U_compare_IF_EB = [];
Q_U_input = [];

%setup for potential step 7 in the future
% equalLines = linspace(-16.25,-7,22);
% comparison = -16.25:0.25:-7;
% diagonalLongitudeList = [];
% 
% for c = 1:length(equalLines)
%     [~,index] = min(abs(comparison-equalLines(c)));
%     diagonalLongitudeList = [diagonalLongitudeList comparison(index)];
% end
% 
% diagonalLatitudeList = 64:-0.25:58.5;
% QU7_collection = [];
% Q7_collection = [];
% Sigma7_collection = [];

for instance = 1:5 %CDO processing creates 5 equally sized files for each variable

   add = num2str(instance);
   fn_U = [u_template add '.nc'];
   
   %may need to multiply by *-1 for some models, CHECK
   U = ncread(fn_U, 'uo');

   Q_V = ncread([Q_template_p1 num2str(instance) Q_template_p2], 'Flow');
   Sigma0 = ncread([sigma0_template_p1 num2str(instance) sigma0_template_p2], 'Sigma0');
   
   if strcmp(modelName, 'CESM2')
       a = (ncread(fn_U, 'time')/365)';
   elseif strcmp(modelName, 'GISS')
       a = (ncread(fn_U, 'time')/365 + 1850)';
   elseif strcmp(modelName, 'IPSL')
       a = (ncread(fn_U, 'time')/365 + 1850)';
   elseif strcmp(modelName, 'EC-Earth3')
       a = (ncread(fn_U, 'time')/365 + 1850)';
   end

   time = [time a]; %check dimensions of time array

   for j = 1:length(a)
        U_specific = squeeze(U(:,:,:,j));
        U_specific = U_specific.*multiplicationMatrix;
        Q_U_gridded = U_specific .* dA./(10^6); %convert to Sv	 
        Q_U_input = cat(4, Q_U_input, Q_U_gridded);%this may need to be changed for parallel computing
        U_gridded = []; 

        Q_V_specific = squeeze(Q_V(:,:,:,j));
        Sigma0_specific = squeeze(Sigma0(:,:,:,j));

        % Step 1: Finding Meridional Transport Across 50N
        Step1 = find(latitude(:,1) == 50); %only focusing on flow on 50N, gulf of st. lawrence corrected
        Q_V_cut = Q_V_specific(Step1,:,:);
        Sigma0_cut = Sigma0_specific(Step1,:,:);
        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_50N = [Flow_50N nansum(nansum(Q_V_cut(mask)))];
        Step1 = [];

        % Step 2: Finding Horizontal Transport at Labrador Sea Western Boundary
        Step2a = find(latitude(:,1) >= 60.5 & latitude(:,1) <= 64); %lat always a due to dimensioning of array
        Step2b = find(longitude(1,:) == -64.5);
        Q_U_gridded_cut = Q_U_gridded(Step2a, Step2b, :);
        Sigma0_cut = Sigma0_specific(Step2a, Step2b, :);
        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_LS_WB = [Flow_LS_WB nansum(nansum(Q_U_gridded_cut(mask)))];
        Step2a = [];
        Step2b = [];

        % Step 3: Finding Meridional Transport at Labrador Sea Northern
        % Boundary
        Step3a = find(latitude(:,1) == 64);
        Step3b = find(longitude(1,:) >= -65 & longitude(1,:) <= -51);
        Q_V_cut = squeeze(Q_V_specific(Step3a, Step3b, :))'; 
        Sigma0_cut = squeeze(Sigma0_specific(Step3a, Step3b, :))';%need to transpose

	    sigma0_compare_LS_NB = cat(3,sigma0_compare_LS_NB, Sigma0_cut);
	    Q_V_compare_LS_NB = cat(3,Q_V_compare_LS_NB, Q_V_cut);
	    long3b = longitude(1,Step3b);

        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_LS_NB = [Flow_LS_NB nansum(nansum(Q_V_cut(mask)))];

        Step3a = [];
        Step3b = [];

        % Step 4: Finding Meridional Transport at Greenland/Iceland Northern Boundary
        Step4a = find(latitude(:,1) == 65.5);
        Step4b = find(longitude(1,:) >= -40 & longitude(1,:) <= -23);
        Q_V_cut = squeeze(Q_V_specific(Step4a, Step4b, :))'; 
        Sigma0_cut = squeeze(Sigma0_specific(Step4a, Step4b, :))';%need to transpose

	    sigma0_compare_GI_NB = cat(3,sigma0_compare_GI_NB, Sigma0_cut);
	    Q_V_compare_GI_NB = cat(3,Q_V_compare_GI_NB, Q_V_cut);
	    long4b = longitude(1,Step4b);

        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_GI_NB = [Flow_GI_NB nansum(nansum(Q_V_cut(mask)))];
        Step4a = [];
        Step4b = [];

        %Steps 5 and 6 might be consolidated in the future
        % Step 5: Finding Meridional Transport at Iceland/Faroe Northern Boundary
        Step5a = find(latitude(:,1) == 64);
        Step5b = find(longitude(1,:) >= -16.5 & longitude(1,:) <= -7);
        Q_V_cut = squeeze(Q_V_specific(Step5a, Step5b, :))'; 
        Sigma0_cut = squeeze(Sigma0_specific(Step5a, Step5b, :))';%need to transpose

	    sigma0_compare_IF_NB = cat(3,sigma0_compare_IF_NB, Sigma0_cut);
	    Q_V_compare_IF_NB = cat(3,Q_V_compare_IF_NB, Q_V_cut);
	    long5b = longitude(1,Step5b);

        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_IF_NB = [Flow_IF_NB nansum(nansum(Q_V_cut(mask)))];
        Step5a = [];
	    Step5b = [];

        % Step 6: Finding Horizontal Transport at Iceland/Faroe Eastern Boundary
        Step6a = find(latitude(:,1) >= 58 & latitude(:,1) <= 64); %lat always 'a'/first due to dimensioning of array
        Step6b = find(longitude(1,:) == -7);
        Q_U_gridded_cut = squeeze(Q_U_gridded(Step6a, Step6b, :))'; %i believe transposing is still necessary here
        Sigma0_cut = squeeze(Sigma0_specific(Step6a, Step6b, :))';

        sigma0_compare_IF_EB = cat(3,sigma0_compare_IF_EB, Sigma0_cut);
	    Q_U_compare_IF_EB = cat(3,Q_U_compare_IF_EB, Q_U_gridded_cut);
	    lat6a = lat(Step6a, 1);

        mask = Sigma0_cut >= meanSigmaMOCvalue;
        Flow_IF_EB = [Flow_IF_EB nansum(nansum(Q_U_gridded_cut(mask)))];
        Step6a = [];
        Step6b = [];

        %Step 7 - Attempted Diagonal
        %for h = 1:length(diagonalLongitudeList)
        %Step7a = find(latitude(:,1) == diagonalLatitudeList(h)); %lat always 'a'/first due to dimensioning of array
        %Step7b = find(long(1,:) == diagonalLongitudeList(h));

        %QU7_collection = [QU7_collection Q_U_gridded(Step7a, Step7b, :)];
        %Q7_collection = [Q7_collection Q_V_specific(Step7a, Step7b, :)];
        %Sigma7_collection = [Sigma7_collection Sigma0_specific(Step7a, Step7b, :)];
        %end
        %mask = Sigma7_collection >= meanSigmaMOCvalue
        %QU7_collection = QU7_collection*-1*sind(45)
        %Q7_collection = Q7_collection*sind(45)*-1
        %Flow_GS_DB = [Flow_GS_DB nansum(nansum(QU7_collection(mask))) + nansum(nansum(Q7_collection(mask)))];
        %mask=[];
        %QU7_collection = [];
        %Q7_collection = [];
	    %Sigma7_collection=[];
   end

%writing .nc file for Q_U for easier future reference
Q_U_filename = ['Q_U_gridded_' modelName '_instance' num2str(instance) '_' num2str(spacing) '_V2.nc'];
[r_nc, c_nc, z_nc, t_nc] = size(Q_U_input);
nccreate(Q_U_filename,'Flow', 'Dimensions',{'x',r_nc,'y', c_nc, 'depth', z_nc, 'time', t_nc});
ncwrite(Q_U_filename, 'Flow', Q_U_input)
Q_U_input = [];

end

%% Format Q_V_compare's and Sigma0_compare's, seasonally filter data of flow
%see helper function below for purpose
sigma0_compare_GI_NB = averageSeasonally(sigma0_compare_GI_NB);
Q_V_compare_GI_NB = averageSeasonally(Q_V_compare_GI_NB);
sigma0_compare_LS_NB = averageSeasonally(sigma0_compare_LS_NB);
Q_V_compare_LS_NB = averageSeasonally(Q_V_compare_LS_NB);
sigma0_compare_IF_NB = averageSeasonally(sigma0_compare_IF_NB);
Q_V_compare_IF_NB = averageSeasonally(Q_V_compare_IF_NB);
sigma0_compare_IF_EB = averageSeasonally(sigma0_compare_IF_EB);
Q_U_compare_IF_EB = averageSeasonally(Q_U_compare_IF_EB);

%might want to rename the above and below functions differently in the
%future for clarity

%% filter inflow data
%trying new annual average, positive is IN FLOW for the northern and
%eastern boundaries
Flow_LS_WB = mean(reshape(Flow_LS_WB, [12 length(Flow_LS_WB)/12]));
Flow_LS_NB = -1*mean(reshape(Flow_LS_NB, [12 length(Flow_LS_NB)/12])); %positive values is net inflow to the control volume
Flow_GI_NB = -1*mean(reshape(Flow_GI_NB, [12 length(Flow_GI_NB)/12]));
Flow_IF_NB = -1*mean(reshape(Flow_IF_NB, [12 length(Flow_IF_NB)/12]));
Flow_IF_EB = -1*mean(reshape(Flow_IF_EB, [12 length(Flow_IF_EB)/12]));
Flow_50N = mean(reshape(Flow_50N, [12 length(Flow_50N)/12])); 

time = floor(min(time)):floor(max(time)); %changing time to annual variable

%% Now time to visualize the IN flow
figure
plot(time, Flow_LS_WB)
hold on
plot(time, Flow_LS_NB) %multiplying by -1 to have consistency
plot(time, Flow_GI_NB)
plot(time, Flow_IF_NB)
plot(time, Flow_IF_EB)
grid on
title(['Annually-Averaged Volumetric INFLOW by Boundary below ' num2str(round(mean(meanSigmaMOCvalue),2))  'kg/m^3 - ' modelName ' (' num2str(spacing) ' spacing)'])
xlabel('Time')
ylabel('Inflow (Sv)')
legend('Lab. Sea. WB.', 'Lab. Sea. NB', 'Greenland/Iceland NB', 'Iceland/Faroe NB', 'Iceland/Faroe EB')
hold off

inputFig1Title = ['VolumeInflow_' modelName '_' num2str(spacing) '_V2.fig'];
savefig(inputFig1Title)

%figure
%plot(time, Flow_LS_WB)
%hold on
%plot(time, Flow_LS_NB*-1) %multiplying by -1 to have consistency
%plot(time, Flow_GI_NB*-1)
%plot(time, Flow_GS_DB)
%grid on
%title(['Volumetric INFLOW by Boundary below ' num2str(round(mean(sigmaMocValue),2))  'kg/m^3 - ' modelName ' (' num2str(spacing) ')'])
%xlabel('Time')
%ylabel('Inflow (Sv)')
%legend('Lab. Sea. WB.', 'Lab. Sea. NB', 'Greenland/Iceland NB', 'Iceland/Faroe DB')
%hold off
%savefig(['VolumeInflow2_compensated_' num2str(spacing) '.fig'])

%a positive totalInflow means water in  the lower limb is flowing into the
%control volume (subpolar north atlantic control area)
totalInflow = Flow_LS_WB + Flow_LS_NB + Flow_GI_NB + Flow_IF_NB + Flow_IF_EB;
Divergence = totalInflow + Flow_50N;

%totalInflow2 = Flow_LS_WB - Flow_LS_NB - Flow_GI_NB + Flow_GS_DB;
%Divergence2 = totalInflow2 + Flow_50N;

figure
plot(time, totalInflow)
hold on
plot(time, Flow_50N)
plot(time, Divergence)
grid on
title(['Annually-Averaged Lower Limb Budget - ' modelName ' (' num2str(spacing) ' spacing)'])
xlabel('Time')
ylabel('Inflow (Sv)')
legend('Total Inflow (positive)', 'Outward Flow at 50N', 'Divergence')
inputFig2Title = ['VolumeBudget_' modelName '_' num2str(spacing) '_V2.fig'];
savefig(inputFig2Title)

input1 = string(['Divergence_Calculations_' modelName '_' num2str(spacing) '_V2.mat']);
input2 = string(['GI_NB_data_' modelName '_' num2str(spacing) '_V2.mat']);
input3 = string(['LS_NB_data_' modelName '_' num2str(spacing) '_V2.mat']);
input4 = string(['IF_NB_data_' modelName '_' num2str(spacing) '_V2.mat']);
input5 = string(['IF_EB_data_' modelName '_' num2str(spacing) '_V2.mat']);

save(input1,getVarName(time),getVarName(Flow_50N),getVarName(Flow_LS_WB),...
    getVarName(Flow_LS_NB),getVarName(Flow_GI_NB),getVarName(Flow_IF_NB), ...
    getVarName(Flow_IF_EB),getVarName(totalInflow),getVarName(Divergence))
save(input2,getVarName(sigma0_compare_GI_NB),getVarName(Q_V_compare_GI_NB),...
    getVarName(long4b),getVarName(depth))
save(input3,getVarName(sigma0_compare_LS_NB),getVarName(Q_V_compare_LS_NB),....
    getVarName(long3b),getVarName(depth))
save(input4,getVarName(sigma0_compare_IF_NB),getVarName(Q_V_compare_IF_NB),...
    getVarName(long5b),getVarName(depth))
save(input5,getVarName(sigma0_compare_IF_EB),getVarName(Q_U_compare_IF_EB),...
    getVarName(lat6a),getVarName(depth))

%% Contour plot for Greenland/Iceland Northern Boundary Section

 input1 = ['GINB_Contour_' modelName '_' num2str(spacing) '_V2.fig'];
 plotMonthlies(Q_V_compare_GI_NB,sigma0_compare_GI_NB, ...
     'Greenland/Iceland Channel Transect @ 65.5N, Averaged Transport vs. Sigma0', ...
     depth, meanSigmaMOCvalue, input1, long4b, 'Longitude (W)')

 input2 = ['LSNB_Contour_' modelName '_' num2str(spacing) '_V2.fig'];
 plotMonthlies(Q_V_compare_LS_NB,sigma0_compare_LS_NB, ...
     'Labrador Sea Northern Transect @ 64N, Averaged Transport vs. Sigma0', ...
     depth, meanSigmaMOCvalue, input2, long3b, 'Longitude (W)')

 input3 = ['IFNB_Contour_' modelName '_' num2str(spacing) '_V2.fig'];
 plotMonthlies(Q_V_compare_IF_NB,sigma0_compare_IF_NB, ...
     'Iceland/Faroe Islands Channel Transect @ 64N, Averaged Transport vs. Sigma0', ...
     depth, meanSigmaMOCvalue, input3, long5b, 'Longitude (W)')

 % still having problems with this
 % input4 = ['IFEB_Contour_' modelName '_' num2str(spacing) '_V2.fig'];
 % plotMonthlies(Q_U_compare_IF_EB,sigma0_compare_IF_NB, ...
 %     'Iceland/Faroe Islands Channel Transect @ -6.75W, Averaged Transport vs. Sigma0', ...
 %     depth, meanSigmaMOCvalue, input4, lat6a, 'Latitude (N)')

function out = getVarName(var)
    out = string(inputname(1));
end

function seasonData_final = averageSeasonally(input)
    %assumes array is originally 3d, [depth, lat/long, time]
     [r,c, t] = size(input);
     monthlyData = reshape(input, [r c 12 t/12]); %only works if time index  is divisible by 
     seasonData_raw = squeeze(mean(monthlyData,4, "omitnan"));

     seasonData_final(:,:,1) = squeeze(mean(seasonData_raw(:,:,1:3), 3)); %January, February, March
     seasonData_final(:,:,2) = squeeze(mean(seasonData_raw(:,:,4:6), 3)); %April, May, June
     seasonData_final(:,:,3) = squeeze(mean(seasonData_raw(:,:,7:9), 3)); %July, August, September
     seasonData_final(:,:,4) = squeeze(mean(seasonData_raw(:,:,10:12), 3)); %October, November, December

     %output array should be a rxcx4 array now
end

function plotMonthlies(Q,Sigma0,FigureTitle, depth, meanSigmaMOCvalue, FigureSave, mapDim, dimFlag)
fig = figure;
subplot(2,2,1)

%find the min sigma in entire dataset
minSigma0 = min([round(min(min(Sigma0(:,:,1), [], 'omitnan')),2) ...
             round(min(min(Sigma0(:,:,2), [], 'omitnan')),2) ...
             round(min(min(Sigma0(:,:,3), [], 'omitnan')),2) ...
             round(min(min(Sigma0(:,:,4), [], 'omitnan')),2)]);

maxSigma0 = max([round(max(max(Sigma0(:,:,1), [], 'omitnan')),2) ...
             round(max(max(Sigma0(:,:,2), [], 'omitnan')),2) ...
             round(max(max(Sigma0(:,:,3), [], 'omitnan')),2) ...
             round(max(max(Sigma0(:,:,4), [], 'omitnan')),2)]);

minFlow = min([round(min(min(Q(:,:,1), [], 'omitnan')),2) ...
             round(min(min(Q(:,:,2), [], 'omitnan')),2) ...
             round(min(min(Q(:,:,3), [], 'omitnan')),2) ...
             round(min(min(Q(:,:,4), [], 'omitnan')),2)]);

minFlow = abs(minFlow); 

maxFlow = max([round(max(max(Q(:,:,1), [], 'omitnan')),2) ...
             round(max(max(Q(:,:,2), [], 'omitnan')),2) ...
             round(max(max(Q(:,:,3), [], 'omitnan')),2) ...
             round(max(max(Q(:,:,4), [], 'omitnan')),2)]);

colorbarAxis = max([minFlow, maxFlow]);

%% find maximum depth where there are still values for each set
depthIndex(1:4) = 0;
[~,c,~] = size(Q);

for z = 1:4
flag = true;
y = 1;
    while flag
        check = sum(isnan(Q(y,:,z)));

        if check == c
        depthIndex(z) = depth(y);
        flag = false;
        end

	y = y+1;
    end
end

depthIndex = min(depthIndex);
depthIndex = find(depth == depthIndex) + 0;

for t = 1:4
    sp{t} = subplot(2,2,t, 'Parent', fig);
    %positive is INFLOW now
    [hh1{t}, hh2(t)] = contourf(sp{t}, mapDim, depth(1:depthIndex), Q(1:depthIndex,:,t)); %need to flip sign
    hh2(t).Color = 'None';
    colormap(othercolor('BuOr_10'))
    caxis(sp{t},[colorbarAxis*-1 colorbarAxis]);

    switch t
        case 1
            subtitle('Jan-Feb-Mar','FontAngle','italic');
        case 2
            subtitle('Apr-May-Jun','FontAngle','italic');
        case 3
            subtitle('Jul-Aug-Sep','FontAngle','italic');
        case 4
            subtitle('Oct-Nov-Dec','FontAngle','italic');
    end

    hold on
    %input = round(linspace(minSigma0+0.01, maxSigma0-0.02,12),2);
    input = [round(minSigma0,1):0.3:maxSigma0-0.05, 27.71];
    %[C{t},h(t)] = contour(mapDim, depth(1:depthIndex), Sigma0(1:depthIndex,:,t), round(linspace(minSigma0+0.05, maxSigma0-0.02,8),2), 'k', "ShowText",true,"LabelFormat","%0.2f");
    %[C{t},h(t)] = contour(mapDim, depth(1:depthIndex), Sigma0(1:depthIndex,:,t), minSigma0+0.05:0.2:maxSigma0-0.02, 'k', "ShowText",true,"LabelFormat","%0.2f");
    [C1{t},h1(t)] = contour(mapDim, depth(1:depthIndex), Sigma0(1:depthIndex,:,t), input(1:2:end-1), 'k');
    [C2{t},h2(t)] = contour(mapDim, depth(1:depthIndex), Sigma0(1:depthIndex,:,t), input(2:2:end), 'k', "ShowText",true,"LabelFormat","%0.2f");
    redLevel = [meanSigmaMOCvalue, meanSigmaMOCvalue];
    [Cr{t},hr(t)] = contour(mapDim,depth(1:depthIndex), Sigma0(1:depthIndex,:,t), redLevel, 'r', 'Linewidth', 1, "ShowText",true,"LabelFormat","%0.2f");
    hold off

end

bv = axes(fig, 'visible', 'off'); 
bv.Title.Visible = 'on';
bv.XLabel.Visible = 'on';
bv.YLabel.Visible = 'on';
ylabel(bv,'Depth (-1*m)','FontSize', 15, 'FontWeight','bold');
xlabel(bv, dimFlag,'FontSize', 15,'FontWeight','bold');
title(bv,FigureTitle,'FontWeight','bold');
c = colorbar(bv,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
caxis(bv,[colorbarAxis*-1 colorbarAxis]);
ylabel(c, 'Inflow Transport [Sv]', 'FontSize', 12,'FontWeight','bold')

savefig(FigureSave)
end
