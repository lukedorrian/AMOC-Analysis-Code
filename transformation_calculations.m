%% Initializations
clear all

addpath GSW %only works if in same higher folder as GSW
addpath GSW/html
addpath GSW/library
addpath GSW/pdf
addpath GSW/thermodynamics_from_t
addpath othercolor
load colorData.mat

%% editable parameters
spacing = 0.5;
modelName = 'GISS';
deltaSigma = 0.2; 

%% more initializations

load_title = string(['regridded_dims_' modelName '_' num2str(spacing) '.mat']);
load(load_title) %loads lat, long saved from calculatingSigmaMOC.m
load_title = string(['sigma_and_MOC_' modelName '_' num2str(spacing) '_ds_' num2str(deltaSigma) '.mat']);
load(load_title)
meanSigmaMOC = round(mean(sigmaMocValue),2);

%Cp = 3850; %J/KgC as per Phys. Oceanography Textbook by UCSD Faculty Lynne Talley
Cp = 4000;
[rGrid,cGrid] = size(long);
tempDx = gsw_distance(long,lat);
tempDx = [tempDx tempDx(:,1)];

dx(1:rGrid,1:cGrid) = tempDx;
dy(1:rGrid,1:cGrid) = gsw_distance([0,0], [0, spacing]);
BF_total = [];
BF_temp_list = [];
BF_Sal_list = [];
BF_total_map = [];
area_list = [];

%loading .nc files
%evaporation and precipitation long/lat variables are not the same, hfds
%long/lat is the same as evaporation (because it's an ocean component)
    if strcmp(modelName, 'CESM2')
        E = ncread('evs_Omon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc', 'evs');
        P = ncread('pr_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','pr'); 
        Q = ncread('hfds_Omon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc', 'hfds');
    elseif strcmp(modelName, 'GISS')
        E = cat(3, ncread('evs_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-190012.nc', 'evs'), ...
            ncread('evs_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_190101-195012.nc', 'evs'), ...
            ncread('evs_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_195101-200012.nc', 'evs'), ...
            ncread('evs_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_200101-201412.nc', 'evs'));

        P = cat(3, ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-190012.nc', 'pr'), ...
            ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_190101-195012.nc', 'pr'), ...
            ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_195101-200012.nc', 'pr'), ...
            ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_200101-201412.nc', 'pr'));

        Q = cat(3, ncread('hfds_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-190012.nc', 'hfds'), ...
            ncread('hfds_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_190101-195012.nc', 'hfds'), ...
            ncread('hfds_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_195101-200012.nc', 'hfds'), ...
            ncread('hfds_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_200101-201412.nc', 'hfds'));

        E = E*-1;
    elseif strcmp(modelName, 'IPSL')
        E = ncread('evs_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc', 'evs');
        P = ncread('pr_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','pr'); 
        Q = ncread('hfds_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc', 'hfds');

        E = E*-1;
    end

Sigma = ncread(['Sigma0_surface_gridded_' modelName '_' num2str(spacing) '.nc'], 'Sigma0');
Theta = ncread(['Theta_surface_gridded_' modelName '_' num2str(spacing) '.nc'], 'thetao');
Sal = ncread(['S_surface_gridded_' modelName '_' num2str(spacing) '.nc'], 'so');

    if strcmp(modelName, 'CESM2')
        latE = ncread('evs_Omon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc', 'lat');
        longE = ncread('evs_Omon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc', 'lon');
        longE(longE > 180) = longE(longE>180) - 360;
        latP = ncread('pr_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','lat');
        longP = ncread('pr_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','lon');
        longP(longP > 180) = longP(longP>180) - 360;
        [longP,latP] = meshgrid(longP,latP);
        longP = longP';
        latP = latP';
    elseif strcmp(modelName, 'GISS')
        %in GISS, coordinates for QEP are all the same
        latP = ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-190012.nc','lat');
        longP = ncread('pr_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-190012.nc','lon');
        longP(longP > 180) = longP(longP>180) - 360;
        [longP,latP] = meshgrid(longP,latP);
        longP = longP';
        latP = latP';
        longE = longP;
        latE = latP;
    elseif strcmp(modelName, 'IPSL')
        latE = double(ncread('evs_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc', 'nav_lat')); %E and Q have same lat/lon
        longE = double(ncread('evs_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc', 'nav_lon'));
        latP = double(ncread('pr_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','lat')); %P is amon so has diff. lat/long from P/Q
        longP = double(ncread('pr_Amon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','lon'));
    end
  

for i = 1:length(time)
E_cut = squeeze(E(:,:,i));
P_cut = squeeze(P(:,:,i));%I really haven't tell if this is accurate...
Q_cut = squeeze(Q(:,:,i));

Sigma_cut = squeeze(Sigma(:,:,i));
Theta_cut = squeeze(Theta(:,:,i));
Sal_cut = squeeze(Sal(:,:,i));

%first, calculate alpha and beta
z_map(1:rGrid,1:cGrid) = 0; %depth array initialization, should I use the first depth?? Ask Yao
press = gsw_p_from_z(z_map,lat); %sea pressure from depth and latitude
sa =  gsw_SA_from_SP(Sal_cut,press,long,lat); %absolute salinity
ct = gsw_CT_from_pt(sa,Theta_cut);
alpha = gsw_alpha(sa,ct,press); %more computationally expensive, maybe more accurate?
beta = gsw_beta(sa,ct,press);

%then regrid E, P, and Q
E_cut = griddata(longE,latE,E_cut, long, lat);
P_cut = griddata(longP,latP,P_cut, long, lat);
Q_cut = griddata(longE,latE,Q_cut, long, lat);

%don't need to mask over P(atmospheric component because final equation
%would apply that mask already
%need mask for control area though
mask = ((lat >= 50 & lat<=65.5) & (long>= -64.5 & long<= -7));
E_cut(~mask) = NaN;
P_cut(~mask) = NaN;
Q_cut(~mask) = NaN;
Sigma_cut(~mask) = NaN;
Theta_cut(~mask) = NaN;
Sal_cut(~mask) = NaN;
alpha(~mask) = NaN;

%two more masks to cut out 64-65.5 latitude areas at Labrador Sea and East
%of Iceland
mask = ((lat > 64 & lat<=65.5) & (long>= -64.5 & long<= -50));
E_cut(mask) = NaN;
P_cut(mask) = NaN;
Q_cut(mask) = NaN;
Sigma_cut(mask) = NaN;
Theta_cut(mask) = NaN;
Sal_cut(mask) = NaN;
alpha(mask) = NaN;
% 
mask = ((lat > 64 & lat<=65.5) & (long>= -16 & long<= -7));
% mask = ((lat > 64 & lat<=65.5) & (long >= -35 & long<= -6.75)); %Only saw
%~10% reduction in heat-induced tranformation
E_cut(mask) = NaN;
P_cut(mask) = NaN;
Q_cut(mask) = NaN;
Sigma_cut(mask) = NaN;
Theta_cut(mask) = NaN;
Sal_cut(mask) = NaN;

%calculate buoyancy transformation, then apply sigma0 mask
sigma0Mask = (meanSigmaMOC - deltaSigma/2 <= Sigma_cut) & (Sigma_cut <= meanSigmaMOC + deltaSigma/2);
BF_temp = dx.*dy.*(-1.*(alpha/Cp).*Q_cut).*(1/10^6).*(1/deltaSigma);
%I made edithere, E_cut + P_cut, I think original equation assumed E_cut
%was positive
BF_Sal = dx.*dy.*beta.*(Sal_cut./(1-Sal_cut)).*(E_cut+P_cut).*(1/10^6).*(1/deltaSigma); %changed E-P to E+P

Area = (dx/1000).*(dy/1000); %converting to km^2
Area(~sigma0Mask) = NaN;
area_list = [area_list round(nansum(nansum(Area)))];

BF_temp(~sigma0Mask) = NaN;
BF_Sal(~sigma0Mask) = NaN;
[temp1,temp2] = size(BF_temp);

BF_temp_list = [BF_temp_list nansum(nansum(BF_temp))];
BF_Sal_list = [BF_Sal_list nansum(nansum(BF_Sal))];
BF_total = [BF_total, nansum(nansum(BF_temp + BF_Sal))];

BF_total_map(1:temp1,1:temp2,i) = BF_temp + BF_Sal;

end

BF_total_map = mean(BF_total_map,3,'omitnan');
%function below assumes vector length is multiple of 12/full years of data
BF_total_annual_avg = mean(reshape(BF_total, [12 length(BF_total)/12]));

input = string(['buoyancyStats_' modelName '_' num2str(spacing) '.mat']);
save(input,"time","BF_temp_list","BF_Sal_list","BF_total", "BF_total_annual_avg", "area_list")

figure
plot(time,BF_total)
grid on
xlabel('Time (years)')
ylabel('Buoyancy Forcing (Sv)')
title(['Buoyancy Forcing over Time - unfiltered-' modelName ', sigmaMoc = ' num2str(meanSigmaMOC) ' kg/m^3, deltaSigma = ' num2str(deltaSigma)])
savefig(['Buoyancy_Transformation_' modelName '_unfiltered_' num2str(spacing) '.fig'])

figure
plot(1849:2013,BF_total_annual_avg)
grid on
xlabel('Time (years)')
ylabel('Annually-Averaged Buoyancy Forcing (Sv)')
title(['Annually-Averaged Buoyancy Forcing over Time - unfiltered-' modelName ', sigmaMoc = ' num2str(meanSigmaMOC) ' kg/m^3, deltaSigma = ' num2str(deltaSigma)])
savefig(['Buoyancy_Transformation_' modelName '_annual_' num2str(spacing) '.fig'])

figure
plot(time,area_list)
grid on
xlabel('Time (years)')
ylabel('Outcropping Area (km^2)')
title('Outcropping Area vs. Time')
savefig(['Area_calculations_' modelName '_' num2str(spacing) '.fig'])

% Past seasonal + 3 year moving filter
% BF_total_seasonal = seasonallyfilterdata(BF_total);
% 
% figure
% plot(time,BF_total_seasonal)
% grid on
% xlabel('Time (years)')
% ylabel('Buoyancy Forcing (Sv)')
% title('Buoyancy Forcing over Time - Seasonal + 3 year moving mean filter')
% savefig(['Buoyancy_Transformation_GISS_filtered_' num2str(spacing) '.fig'])

figure
contourf(long,lat,BF_total_map)
colormap(othercolor('BuOr_10'))
hold on
Sigma_cut_copy = Sigma_cut;
Sigma_cut_copy(~isnan(Sigma_cut_copy)) = max(max(BF_total_map))+0.01;
Sigma_cut_copy(isnan(Sigma_cut_copy)) = 0;
b = colorbar;
ylabel(b,'Transformation (Sv)','FontSize',10,'Rotation',270);
contour(long,lat,Sigma_cut_copy, [0,max(max(BF_total_map))+0.01], 'k', 'Linewidth', 1.5)
title(['Averaged Surface Buoyancy Transformation, ' modelName ', ' num2str(spacing)  ' Spacing, Max Trans. Rate: ' num2str(round(max(max(BF_total_map)))) ' Sv'])
xlabel('Longitude')
ylabel('Latitude')
inputFig1Title = ['ABT_' modelName '_' num2str(spacing) '.fig'];
savefig(inputFig1Title)
hold off

% figure
% contourf(long,lat,Sigma_cut)
% colormap(othercolor('BuOr_10'))
% hold on
% Sigma_cut_copy = Sigma_cut;
% Sigma_cut_copy(~isnan(Sigma_cut_copy)) = max(max(Sigma_cut))+0.1;
% Sigma_cut_copy(isnan(Sigma_cut_copy)) = min(min(Sigma_cut))-0.1;
% a = colorbar;
% ylabel(a,'Sigma0 (kg/m^3)','FontSize',10,'Rotation',270);
% contour(long,lat,Sigma_cut_copy, [min(min(Sigma_cut))-0.1,max(max(Sigma_cut))+0.1;], 'k', 'Linewidth', 1.5)
% contourf(long,lat,Sigma_cut, [meanSigmaMOC - deltaSigma/2,meanSigmaMOC + deltaSigma/2;], 'g--')
% title(['Sigma0 Surface Values, ' modelName ', ' num2str(spacing)  ' Spacing, January 1850'])
% xlabel('Longitude')
% ylabel('Latitude')
% inputFig2Title = ['Sigma_Surface_' modelName '_' num2str(spacing) '.fig'];
% savefig(inputFig2Title)
% hold off

% function filtered = seasonallyfilterdata(Data) %assuming a row vector or a 2D vector
% 
%     for i = 1:12
%     Data(:,i:12:end) = Data(:,i:12:end) - mean(Data(:,i:12:end)) + mean(Data);
%     end
%     filtered = movmean(Data,3*12);
% end

function out = getVarName(var)
    out = string(inputname(1));
end