% What does this code do?
% 1. Computes meridional volumetric transfer from model velocities 
% 2. Computes Sigma0 from Salinity and Theta
% 3. Saves regridded meridional Q velocities as an .nc file
% 4. Saves regridded surface Salinity and Theta values for quicker
% transformation calculations
% 5. Plots meridional volumetric transfer across 50 N
% 6. Calculates and plots Sigma_MOC (Sigma0 value that differentiates upper
% and lower limb) over time, calculates average based on density bins 
% Note: Sigma_MOC value differs based on density bin established
% You will need to rerun this code to get data for different spacing,or
% densitySpacing values. This code provides other data that is calculated
% with divergence_calculations.m and transformation_calculations.m
clear all

%%initializations
addpath GSW %only works if in same higher folder as GSW
addpath GSW/html
addpath GSW/library
addpath GSW/pdf
addpath GSW/thermodynamics_from_t

%% ADJUSTABLE PARAMETERS, EDIT HERE TO CHANGE MODEL, SPACING, AND DELTA SIGMA
modelName = 'CESM2'; 
spacing = 0.5; %in degrees longitude/latitude
densitySpacing = 0.2; %in kg/m^3 - using 0.2 kg/m^3 because this because Speer et. al. 1992 claims it will reduce noise at finer intervals

%% relevant .nc files and other parameters

if strcmp(modelName, 'GISS')
    s_template = 'so_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_'; %need to add number + .nc'
    theta_template = 'thetao_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_'; %need to add number + .nc'
    v_template = 'vo_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_'; %need to add number + .nc'
    s_example = 'so_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-187012.nc';
    v_example = 'vo_Omon_GISS-E2-1-G_historical_r1i1p1f1_gn_185001-187012.nc';
elseif strcmp(modelName, 'CESM2')
    s_template = 'so_Omon_CESM2_historical_r1i1p1f1_gn_ATL00000'; %need to add number + .nc'
    theta_template = 'thetao_Omon_CESM2_historical_r1i1p1f1_gn_ATL00000'; %need to add number + .nc'
    v_template = 'vo_Omon_CESM2_historical_r1i1p1f1_gn_ATL00000'; %need to add number + .nc'
    s_example = 'so_Omon_CESM2_historical_r1i1p1f1_gn_ATL000000.nc';
    v_example = 'vo_Omon_CESM2_historical_r1i1p1f1_gn_ATL000000.nc';
elseif strcmp(modelName, 'IPSL')
    %so_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-194912.nc
    %so_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_195001-201412.nc
    s_template = 'so_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_';
    theta_template = 'thetao_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_'; %need to add number + .nc'
    v_template = 'vo_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_'; %need to add number + .nc'
    s_example = 'so_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-194912.nc';
    v_example = 'vo_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-194912.nc';

end

%% retrieve and adjust dimensional coordinates
%longitude and latitude coordinates for salinity/theta data calculations
%these are usually model-specific!

%making mask since only currently concerned about region around 50N
%maybe just regrid all data longitude -55 to -5
%latitude, 49 to 52
longitudeRange = -65:spacing:-5;
latitudeRange = 49:spacing:67;
[long,lat] = meshgrid(longitudeRange,latitudeRange);
R = find(lat(:,1) == 50); %only focusing on flow on 50N

%adding lat and long mask for small Gulf of St. Lawrence Correction
latMask = lat >= 49 & lat <= 52.25;
longMask = long>= -65 & long <= -56;
combinedMask = latMask & longMask;
latMask = [];
longMask = [];

%longitude and latitude for velocity dataset, model-specific
if strcmp(modelName, 'CESM2')||strcmp(modelName, 'GISS')
longitude = ncread(v_example, 'lon'); %convert to degrees east, what is longitude variable in .nc file?
longitude(longitude > 180) = longitude(longitude>180) - 360;
latitude = ncread(v_example, 'lat'); %degrees north, what is latitude variable in .nc file?

longS = ncread(s_example, 'lon'); %longS doesn't match long of model sometimes
longS(longS > 180) = longS(longS>180) - 360;
latS = ncread(s_example, 'lat'); %matches theta

    %GISS-specific adjustments to associated variables in this if statement
    if strcmp(modelName, 'GISS')
    [longS,latS] = meshgrid(longS,latS);
    longS = longS';
    latS = latS';
    [rS,cS] = size(latS);
    %more GISS-specific changes
    [longitude,latitude] = meshgrid(longitude,latitude);
    longitude = longitude';
    latitude = latitude';
    end

elseif strcmp(modelName, 'IPSL')
longitude = double(ncread(v_example, 'nav_lon')); %convert to degrees east, what is longitude variable in .nc file?
%longitude(longitude > 180) = longitude(longitude>180) - 360; Don't think I
%need to do this conversion, might want to double check
latitude = double(ncread(v_example, 'nav_lat')); %degrees north, what is latitude variable in .nc file?

longS = double(ncread(s_example, 'nav_lon')); %longS doesn't match long of model sometimes
%longS(longS > 180) = longS(longS>180) - 360;
latS = double(ncread(s_example, 'nav_lat')); %matches theta
end

[rS,cS] = size(latS);

%depth-specific reading changes, model-specific
if strcmp(modelName, 'GISS')
    depth = ncread(v_example, 'lev')/(-1); %changed for GISS
    d_bounds = ncread(v_example, 'lev_bnds');
elseif strcmp(modelName, 'CESM2')
    depth = ncread(v_example, 'lev')/(-100); %Convert from cm to meters
    d_bounds = ncread(v_example, 'lev_bnds');
elseif strcmp(modelName, 'IPSL')
    depth = double(ncread(v_example, 'olevel'))/(-1);
    d_bounds = double(ncread(v_example, 'olevel_bounds'));
end

D = length(depth);
delta_z = (d_bounds(2,:) - d_bounds(1,:))'; %STILL APPLIES TO GISS,CESM2,IPSL but might change for other models

v_example = [];
s_example = [];

[rGrid,cGrid] = size(long);
dA(1:rGrid,1:cGrid,1:D) = 0;

for z = 1:length(depth)
    z_map(1:rGrid,1:cGrid) = depth(z); %depth array initialization
    p = gsw_p_from_z(z_map,lat); %sea pressure from depth and latitude
    tempDx = gsw_distance(long,lat,p);
    tempDx = [tempDx tempDx(:,end)];
    dx(1:rGrid,1:cGrid) = tempDx; 
    dA(:,:,z) = dx*delta_z(z);
end

%% initializing variables
z_map = [];
add = [];
time = [];
Q_density = []; 
sigmaMocValue = [];
max_Q_value = [];

Q_input_nc = [];
sigma0_input_nc = [];
sigma0_surface_nc = [];
sal_surface_nc = [];
theta_surface_nc = [];

densityRange = 24:densitySpacing:29; %potential density range

p(1:rS,1:cS) = 0; %initialize pressure array/ No reason to store this data as a 3D array
sa = p; %initialize absolute salinity array (for in-situ temperature calculation)
ct = p; %initialize in situ temperature array (calculated from potential temperature array)
%^change name if necessary

%determine how many iterations are needed based on number of model .nc
%files for each relevant variable
if strcmp(modelName, 'CESM2')
lastFileNum = 3;
elseif strcmp(modelName, 'GISS')
lastFileNum = 8;
elseif strcmp(modelName, 'IPSL')
lastFileNum = 1;
end

%% beginning of loop
for instance = 0:lastFileNum
    switch instance
        case 0
            if strcmp(modelName, 'GISS')
                add = '185001-187012';
            elseif strcmp(modelName, 'CESM2')
                add = '0';
            elseif strcmp(modelName, 'IPSL')
                add = '185001-194912';
            end
        case 1
            if strcmp(modelName, 'GISS')
                add = '187101-189012';
            elseif strcmp(modelName, 'CESM2')
                add = '1';
            elseif strcmp(modelName, 'IPSL')
                add = '195001-201412';
            end
        case 2
            if strcmp(modelName, 'GISS')
                add = '189101-191012';
            elseif strcmp(modelName, 'CESM2')
                add = '2';
            end
        case 3
            if strcmp(modelName, 'GISS')
                add = '191101-193012';
            elseif strcmp(modelName, 'CESM2')
                add = '3';
            end
        case 4
            add = '193101-195012';
        case 5
            add = '195101-197012';
        case 6 
            add = '197101-199012';
        case 7
            add = '199101-201012';
        case 8 
            add = '201101-201412';
    end

   fn_V = [v_template add '.nc'];
   fn_sal = [s_template add '.nc'];
   fn_theta = [theta_template add '.nc'];
   
   V = ncread(fn_V, 'vo');
   S = ncread(fn_sal, 'so');
   theta = ncread(fn_theta, 'thetao');

   if strcmp(modelName, 'GISS')
       a = (ncread(fn_V, 'time')/365 + 1850)';
   elseif strcmp(modelName, 'CESM2')
       a = (ncread(fn_V, 'time')/365)';
   elseif strcmp(modelName, 'IPSL')
       V = double(V); %IPSL-specific corrections since model data is a single
       S = double(S);
       theta=double(theta);
       a = (ncread(fn_V, 'time')/365 + 1850)';
   end

%reducing computing time, taking every third entry
   time = [time a];

   for i = 1:length(a) %reducing computing time
        S_specific = squeeze(S(:,:,:,i));
        theta_specific = squeeze(theta(:,:,:,i));
        [rt,ct,zt]= size(theta_specific);
        sigma0_specific(1:rt,1:ct,1:zt) = 0;
        rt= [];
        ct = [];
        zt = [];

            for z = 1:D
                  z_map(1:rS,1:cS) = depth(z); %depth array initialization
                  p = gsw_p_from_z(z_map,latS); %sea pressure from depth and latitude
                  sa =  gsw_SA_from_SP(S_specific(:,:,z),p,longS,latS); %absolute salinity
                  ct = gsw_CT_from_pt(sa,theta_specific(:,:,z));
                  sigma0_specific(:,:,z) = gsw_sigma0(sa,ct);
            end
        
        V_specific = squeeze(V(:,:,:,i));%loading later to save memory
        V_gridded(1:rGrid,1:cGrid,1:length(depth)) = 0;
        sigma0_gridded = V_gridded;

         for z = 1:D
              V_gridded(:,:,z) = griddata(longitude,latitude,V_specific(:,:,z), long, lat);
              sigma0_gridded(:,:,z) = griddata(longS,latS,sigma0_specific(:,:,z), long, lat);
                
              %section added to save regridded surface and theta values
              %for transformation calculations
              if z == 1
              S_gridded = griddata(longS,latS,S_specific(:,:,z), long,lat);
              Theta_gridded = griddata(longS,latS,theta_specific(:,:,z), long,lat);

              sal_surface_nc = cat(3,sal_surface_nc, S_gridded);
              theta_surface_nc = cat(3,theta_surface_nc, Theta_gridded);
              end

		      %adding small correction for Gulf of St. Lawrence
	          V_copy = V_gridded(:,:,z);
	          V_copy(combinedMask) = NaN;
	          V_gridded(:,:,z) = V_copy;
	          V_copy = [];
         end
         
         S_specific = [];
         theta_specific = [];
         V_specific = [];
         Q = V_gridded.*dA./(10^6); %calculate flow for entire dataset, in Sv (10^6/s)
         V_gridded = [];

	%compensation method below
	[rQ,~,~] = size(Q);

	for f = 1:rQ
	
	    %I am finding the discrepancy of the cumulative volume transfer at the bottom, then compensating by adding flow/number of cells with non-NaN values

	    Q_cut = squeeze(Q(f,:,:));
	    Q_cut_copy = Q_cut; %copy for figuring out how many non NaN-s to divide among 
	    maskNaN = isnan(Q_cut_copy);
	    Q_cut_copy(~isnan(Q_cut_copy)) = 1;
	    Q_cut_copy(isnan(Q_cut_copy)) = 0;
	    sumReal = sum(sum(Q_cut_copy));

	    Q_cut_copy = Q_cut;
	    Q_cut_copy = nansum(Q_cut_copy,1);%nansum(array,2) was WRONG. Will probably need to run through calculcations again......
	    Q_cut_copy = cumsum(Q_cut_copy);
	    difference = Q_cut_copy(end); %Theoretically, this should be 0 but usually isn't due to resolution issues and average flows within a specific cell
	    offset = difference/sumReal;
	    Q_cut(~isnan(Q_cut)) = Q_cut(~isnan(Q_cut)) - offset; %offsets all non-NaN volumetric transfer cells along a certain latitude
	    Q(f,:,:) = Q_cut; %corrects the original volumetric flow dataset by offset 
	end

	Q_input_nc = cat(4,Q_input_nc,Q);
        Q_cut = squeeze(Q(R,:,:)); %had a 3d array around 50N, now summing across just 50N
        Q = [];
        sigma0_cut = squeeze(sigma0_gridded(R,:,:));
	sigma0_input_nc = cat(4, sigma0_input_nc,sigma0_gridded);
         
    %adding part for sigma surface .nc file
    sigma0_surface_nc = cat(3, sigma0_surface_nc,squeeze(sigma0_gridded(:,:,1)));
    
          for x = 1:length(densityRange)
                mask = densityRange(x) - densitySpacing/2 <= sigma0_cut & sigma0_cut <= densityRange(x) + densitySpacing/2; 
                Qvalue = nansum(Q_cut(mask));
                Q_density(x) = Qvalue;
          end
            
          fi_density = cumsum(Q_density); %easy because this is a vector now
          [max_Q_value_1,index] = max(fi_density);
          max_Q_value = [max_Q_value max_Q_value_1];
          sigmaMocValue = [sigmaMocValue, densityRange(index)];
         
   end

   V = [];
   S = [];
   theta = [];

[r_nc,c_nc,z_nc,t_nc] = size(Q_input_nc);

%first write gridded Q .nc file
Q_filename = ['Q_V_gridded_' modelName '_instance' num2str(instance) '_' num2str(spacing) '.nc'];
disp(Q_filename)
nccreate(string(Q_filename),"Flow", ...
"Dimensions",{"x",r_nc,"y", c_nc, "depth", z_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(Q_filename), "Flow", Q_input_nc)

Q_input_nc = [];
Q_filename = [];

sigma0_filename = ['Sigma0_gridded_' modelName '_instance' num2str(instance) '_' num2str(spacing) '.nc'];
nccreate(string(sigma0_filename), "Sigma0", ...
"Dimensions", {"x", r_nc, "y", c_nc, "depth", z_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(sigma0_filename), "Sigma0", sigma0_input_nc)

sigma0_input_nc = [];
sigma0_filename = [];
%keeping .nc file extremely simple right now as I had previous problems with nccreate(), going to save other important details about the regridded file as a .mat file currently
end

%write .nc file for surface sigma0, sal, theta values
[r_nc,c_nc,t_nc] = size(sigma0_surface_nc);
sigma0_filename = ['Sigma0_surface_gridded_' modelName '_' num2str(spacing) '.nc'];
nccreate(string(sigma0_filename), "Sigma0", ...
"Dimensions", {"x", r_nc, "y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(sigma0_filename), "Sigma0", sigma0_surface_nc)

S_filename = ['S_surface_gridded_' modelName '_' num2str(spacing) '.nc'];
nccreate(string(S_filename),"so", "Dimensions",{"x",r_nc,"y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(S_filename), "so", sal_surface_nc)

Theta_filename = ['Theta_surface_gridded_' modelName '_' num2str(spacing) '.nc'];
nccreate(string(Theta_filename),"thetao", "Dimensions",{"x",r_nc,"y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(Theta_filename), "thetao", theta_surface_nc)

%% Calculating Annual Average MOC Strength and Average Sigma MOC
annual_max_Q_value = mean(reshape(max_Q_value, [12 length(max_Q_value)/12]));
meanSigmaMoc = mean(sigmaMocValue);

%% Figures

figure
plot(time, sigmaMocValue, 'b', 'Linewidth', 1.5)
hold on
plot(time, meanSigmaMoc, 'r')
hold off
grid on
xlabel('Time [years]')
ylabel('Potential Density (kg/m^3)')
title(['Sigma MOC over Time - Model: '  modelName])
legend('SigmaMOC(t)', 'Avg. SigmaMOC')
savefig(['SigmaMoc_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '.fig'])

figure
plot(time, max_Q_value)
grid on
xlabel('Time [years]')
ylabel('Maximum Cumulative Flow (Sv) ')
title(['Maximum Cumulative Flow over 50N - Model: '  modelName ', Spacing: ' num2str(spacing)])
savefig(['max_Qflow_50N_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '.fig'])

figure
plot(1849:2013, annual_max_Q_value, 'b', 'Linewidth', 1.5)
grid on
xlabel('Time [years]')
ylabel('Maximum Cumulative Flow (Sv) ')
title(['Annual Average Maximum Cumulative Flow over 50N - Model: '  modelName ', Spacing: ' num2str(spacing)])
savefig(['annual_avg_Qflow_50N_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '.fig'])

%make sure to rename these when starting a new model!
input1 = string(['sigma_and_MOC_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '.mat']);
input2 = string(['regridded_dims_' modelName '_' num2str(spacing) '.mat']);

save(input1, getVarName(time), getVarName(max_Q_value), getVarName(sigmaMocValue), ...
    getVarName(densitySpacing), getVarName(annual_max_Q_value)) 
save(input2, getVarName(time), getVarName(lat), getVarName(long), ...
    getVarName(depth), getVarName(delta_z))

function out = getVarName(var)
    out = string(inputname(1));
end
