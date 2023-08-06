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
modelName = 'CESM2';
deltaSigma = 0.2; 
%% more initializations

load_title = string(['regridded_dims_' modelName '_' num2str(spacing) '_V2.mat']);
load(load_title) %loads latitude, longitude saved from calculatitudeingSigmaMOC.m
load_title = string(['sigma_and_MOC_' modelName '_' num2str(spacing) '_ds_' num2str(deltaSigma) '_V2.mat']);
load(load_title)
meanSigmaMOC = round(mean(sigmaMocValue),2);

Cp = 4000; %specific heat capacity of seawater
[rGrid,cGrid] = size(longitude);
tempDx = gsw_distance(longitude,latitude);
tempDx = [tempDx tempDx(:,1)];

dx(1:rGrid,1:cGrid) = tempDx;
dy(1:rGrid,1:cGrid) = gsw_distance([0,0], [0, spacing]);
BF_total = [];
BF_temp_list = [];
BF_Sal_list = [];
BF_total_map = [];
area_list = [];

%% loading .nc files
%evaporation and precipitation longitude/latitude variables are not the same, hfds
%longitude/latitude is the same as evaporation (because it's an ocean component)

%now, all files are already regridded via CDO and  combined. They have the
%same longitude/latitudeitude coordinates and the same name templatitudee

E = ncread(['evs_' modelName '_regridded.nc'], 'evs');
P = ncread(['pr_' modelName '_regridded.nc'],'pr'); 
Q = ncread(['hfds_' modelName '_regridded.nc'], 'hfds');

%!!! - Model-Specific Change
if strcmp(modelName,'GISS')||strcmp(modelName,'IPSL')||strcmp(modelName,'EC-Earth3')
E = E*-1;
end

Sigma = ncread(['Sigma0_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'], 'Sigma0');
Theta = ncread(['Theta_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'], 'thetao');
Sal = ncread(['S_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'], 'so');

%the latitudeitudes and longitudes are all the same now because all data has now
%been regridded
  
%% loop part

for i = 1:length(time)
E_cut = squeeze(E(:,:,i));
P_cut = squeeze(P(:,:,i));
Q_cut = squeeze(Q(:,:,i));

Sigma_cut = squeeze(Sigma(:,:,i));
Theta_cut = squeeze(Theta(:,:,i));
Sal_cut = squeeze(Sal(:,:,i));

%first, calculate alpha and beta
z_map(1:rGrid,1:cGrid) = 0; %depth array initialization, should I use the first depth?? Ask Yao
press = gsw_p_from_z(z_map,latitude); %sea pressure from depth and latitudeitude
sa =  gsw_SA_from_SP(Sal_cut,press,longitude,latitude); %absolute salinity
ct = gsw_CT_from_pt(sa,Theta_cut);
alpha = gsw_alpha(sa,ct,press); %more computationally expensive, maybe more accurate?
beta = gsw_beta(sa,ct,press);

%need mask for control area though
mask = ((latitude >= 50 & latitude<=65.5) & (longitude>= -64.5 & longitude<= -7));
E_cut(~mask) = NaN;
P_cut(~mask) = NaN;longitude
Q_cut(~mask) = NaN;
Sigma_cut(~mask) = NaN;
Theta_cut(~mask) = NaN;
Sal_cut(~mask) = NaN;
alpha(~mask) = NaN;

%two more masks to cut out 64-65.5 latitudeitude areas at Labrador Sea and East
%of Iceland
mask = ((latitude > 64 & latitude<=65.5) & (longitude>= -64.5 & longitude<= -50));
E_cut(mask) = NaN;
P_cut(mask) = NaN;
Q_cut(mask) = NaN;
Sigma_cut(mask) = NaN;
Theta_cut(mask) = NaN;
Sal_cut(mask) = NaN;
alpha(mask) = NaN;
% 
mask = ((latitude > 64 & latitude<=65.5) & (longitude>= -16 & longitude<= -7));

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
%line below assumes vector length is multiple of 12/full years of data
BF_total_annual_avg = mean(reshape(BF_total, [12 length(BF_total)/12]));

input = string(['buoyancyStats_' modelName '_' num2str(spacing) '_V2.mat']);
save(input,"time","BF_temp_list","BF_Sal_list","BF_total", "BF_total_annual_avg", "area_list")

figure
plot(time,BF_total)
grid on
xlabel('Time (years)')
ylabel('Buoyancy Forcing (Sv)')
title(['Buoyancy Forcing over Time - unfiltered-' modelName ', sigmaMoc = ' num2str(meanSigmaMOC) ' kg/m^3, deltaSigma = ' num2str(deltaSigma)])
savefig(['Buoyancy_Transformation_' modelName '_unfiltered_' num2str(spacing) '_V2.fig'])

figure
plot(1850:2014,BF_total_annual_avg)
grid on
xlabel('Time (years)')
ylabel('Annually-Averaged Buoyancy Forcing (Sv)')
title(['Annually-Averaged Buoyancy Forcing over Time - unfiltered-' modelName ', sigmaMoc = ' num2str(meanSigmaMOC) ' kg/m^3, deltaSigma = ' num2str(deltaSigma)])
savefig(['Buoyancy_Transformation_' modelName '_annual_' num2str(spacing) '_V2.fig'])

figure
plot(time,area_list)
grid on
xlabel('Time (years)')
ylabel('Outcropping Area (km^2)')
title('Outcropping Area vs. Time')
savefig(['Area_calculatitudeions_' modelName '_' num2str(spacing) '_V2.fig'])

figure
contourf(longitude,latitude,BF_total_map)
colormap(othercolor('BuOr_10'))
hold on
Sigma_cut_copy = Sigma_cut;
Sigma_cut_copy(~isnan(Sigma_cut_copy)) = max(max(BF_total_map))+0.01;
Sigma_cut_copy(isnan(Sigma_cut_copy)) = 0;
b = colorbar;
ylabel(b,'Transformation (Sv)','FontSize',10,'Rotation',270);
contour(longitude,latitude,Sigma_cut_copy, [0,max(max(BF_total_map))+0.01], 'k', 'Linewidth', 1.5)
title(['Averaged Surface Buoyancy Transformation, ' modelName ', ' num2str(spacing)  ' Spacing, Max Trans. Rate: ' num2str(round(max(max(BF_total_map)))) ' Sv'])
xlabel('longitude')
ylabel('latitudeitude')
inputFig1Title = ['ABT_' modelName '_' num2str(spacing) '_V2.fig'];
savefig(inputFig1Title)
hold off

% figure
% contourf(longitude,latitude,Sigma_cut)
% colormap(othercolor('BuOr_10'))
% hold on
% Sigma_cut_copy = Sigma_cut;
% Sigma_cut_copy(~isnan(Sigma_cut_copy)) = max(max(Sigma_cut))+0.1;
% Sigma_cut_copy(isnan(Sigma_cut_copy)) = min(min(Sigma_cut))-0.1;
% a = colorbar;
% ylabel(a,'Sigma0 (kg/m^3)','FontSize',10,'Rotation',270);
% contour(longitude,latitude,Sigma_cut_copy, [min(min(Sigma_cut))-0.1,max(max(Sigma_cut))+0.1;], 'k', 'Linewidth', 1.5)
% contourf(longitude,latitude,Sigma_cut, [meanSigmaMOC - deltaSigma/2,meanSigmaMOC + deltaSigma/2;], 'g--')
% title(['Sigma0 Surface Values, ' modelName ', ' num2str(spacing)  ' Spacing, January 1850'])
% xlabel('longitude')
% ylabel('latitudeitude')
% inputFig2Title = ['Sigma_Surface_' modelName '_' num2str(spacing) '.fig'];
% savefig(inputFig2Title)
% hold off

function out = getVarName(var)
    out = string(inputname(1));
end