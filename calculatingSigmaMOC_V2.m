%Rewriting Code in lieu of new CDO method that unifies input data stream
% What does this code do?
% 1. Computes meridional volumetric transfer from model velocities 
% 2. Computes Sigma0 from Salinity and Theta
% 3. Saves regridded meridional Q velocities as an .nc file
% 4. Saves regridded surface Salinity and Theta values for quicker transformation calculations
% 5. Plots meridional volumetric transfer across 50 N
% 6. Calculates and plots Sigma_MOC (Sigma0 value that differentiates upper
% and lower limb) over time, calculates average based on density bins 

% Note: Sigma_MOC value differs based on density bin established
% You will need to rerun this code to get data for different spacing,or
% densitySpacing values. This code provides other data that is calculated
% with divergence_calculations.m and transformation_calculations.m

%% initializations
%clear all
addpath GSW %only works if in same higher folder as GSW
addpath GSW/html
addpath GSW/library
addpath GSW/pdf
addpath GSW/thermodynamics_from_t

%% ADJUSTABLE PARAMETERS, EDIT HERE TO CHANGE MODEL, SPACING, AND DELTA SIGMA
modelName = 'CESM2'; 
spacing = 0.5; %MAKE SURE TO ENSURE THIS MATCHES GRID OF INPUT DATA, This does not change the input data, just adjusts figure titles
densitySpacing = 0.2; %in kg/m^3 - using 0.2 kg/m^3 because this because Speer et. al. 1992 claims it will reduce noise at finer intervals
densityRange = 24:densitySpacing:29; %potential density range

%% relevant .nc files and other parameters
s_template = ['so_' modelName '_regridded_subset_'];
theta_template = ['thetao_' modelName '_regridded_subset_'];
v_template = ['vo_' modelName '_regridded_subset_'];
s_example = ['so_' modelName '_regridded_subset_000001.nc'];
v_example = ['vo_' modelName '_regridded_subset_000001.nc'];

%% Reading latitude,longitude and depth coordinates, non-model specific (with CDO process)
%CDO formatting appears to change the name of the latitude and longitude
%coordinates to make them all the same, now since regridding has already
%occurred, S,Theta coordinates match with U,V

%for details on the gridding, see grid_ATL_50km.txt (for 0.5 spacing)
longitude = ncread(v_example, 'LONGITUDE'); %convert to degrees east, what is longitude variable in .nc file?
latitude = ncread(v_example, 'LATITUDE'); %degrees north, what is latitude variable in .nc file?

[longitude,latitude] = meshgrid(longitude,latitude);
longitude = longitude';
latitude = latitude';

%adding lat and long mask for small Gulf of St. Lawrence Correction
latMask = latitude >= 49 & latitude <= 52.25;
longMask = longitude >= -65 & longitude <= -56;
combinedMask = latMask & longMask;
latMask = [];
longMask = [];
R = find(latitude(1,:) == 50); %only focusing on calculating flow on 50N for this code

%% !!!!! MODEL-SPECIFIC changes for depth, need to write code here when applying to new model set
%CDO does NOT change depth variable name and values based on current method
if strcmp(modelName, 'GISS')||strcmp(modelName,'EC-Earth3')
    depth = ncread(v_example, 'lev')/(-1); %changed for GISS
    d_bounds = ncread(v_example, 'lev_bnds');
elseif strcmp(modelName, 'CESM2')
    depth = ncread(v_example, 'lev')/(-100); %Convert from cm to meters
    d_bounds = ncread(v_example, 'lev_bnds');
elseif strcmp(modelName, 'IPSL')
    depth = double(ncread(v_example, 'olevel'))/(-1);
    d_bounds = double(ncread(v_example, 'olevel_bounds'));
end

%% non-model-specific code
D = length(depth);

delta_z = (d_bounds(2,:) - d_bounds(1,:))'; %STILL APPLIES TO GISS,CESM2,IPSL but might change for other models

[rGrid,cGrid] = size(longitude);
combinedMask = repmat(combinedMask, [1, 1, D]);
multiplicationMatrix(1:rGrid,1:cGrid,1:length(depth)) = 1;
multiplicationMatrix(combinedMask)=NaN;
dA(1:rGrid,1:cGrid,1:D) = 0;

for z = 1:length(depth)
    z_map(1:rGrid,1:cGrid) = depth(z); %depth array initialization
    p = gsw_p_from_z(z_map,latitude); %sea pressure from depth and latitude
    tempDx = gsw_distance(longitude,latitude,p);
    tempDx = [tempDx tempDx(:,end)];
    dx(1:rGrid,1:cGrid) = tempDx; 
    dA(:,:,z) = dx*delta_z(z);
end

%% Setting up variables for the loop
v_example = [];
s_example = [];
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

p(1:rGrid,1:cGrid) = 0; %initialize pressure array/ No reason to store this data as a 3D array
SA = p; %initialize absolute salinity array (for in-situ temperature calculation)
CT = p; %initialize in situ temperature array (calculated from potential temperature array)

%% beginning of loop
for instance = 1:5

    add = ['00000' num2str(instance)];

   fn_V = [v_template add '.nc'];
   fn_sal = [s_template add '.nc'];
   fn_theta = [theta_template add '.nc'];

   %the class name for each variable might be model-specific,may need to adjust in the future
   V = ncread(fn_V, 'vo');
   S = ncread(fn_sal, 'so');
   theta = ncread(fn_theta, 'thetao');

   %model-specific time changes
   if strcmp(modelName, 'GISS')
       a = (ncread(fn_V, 'time')/365 + 1850)';
   elseif strcmp(modelName, 'CESM2')
       a = (ncread(fn_V, 'time')/365)';
   elseif strcmp(modelName, 'IPSL')
       V = double(V); %IPSL-specific corrections since model data is a single
       S = double(S);
       theta=double(theta);
       a = (ncread(fn_V, 'time')/365 + 1850)';
   elseif strcmp(modelName,'EC-Earth3')
       a = (ncread(fn_V, 'time')/365 + 1850)';
   end

   time = [time a];

   for i = 1:length(a)
        S_specific = squeeze(S(:,:,:,i));
        theta_specific = squeeze(theta(:,:,:,i));
        V_specific = squeeze(V(:,:,:,i));
        V_specific = V_specific.*multiplicationMatrix; %gulf of st.lawrence correction
        sigma0_specific(1:rGrid,1:cGrid,1:D) = 0;

            for z = 1:D
                  z_map(1:rGrid,1:cGrid) = depth(z); %depth array initialization
                  p = gsw_p_from_z(z_map,latitude); %sea pressure from depth and latitude
                  SA =  gsw_SA_from_SP(S_specific(:,:,z),p,longitude,latitude); %absolute salinity
                  CT = gsw_CT_from_pt(SA,theta_specific(:,:,z));
                  sigma0_specific(:,:,z) = gsw_sigma0(SA,CT);
            end

        sal_surface_nc = cat(3,sal_surface_nc, squeeze(S_specific(:,:,1)));
        theta_surface_nc = cat(3,theta_surface_nc, squeeze(theta_specific(:,:,1)));
        sigma0_surface_nc = cat(3, sigma0_surface_nc,squeeze(sigma0_specific(:,:,1)));

         Q = V_specific.*dA./(10^6); %calculate flow for entire dataset, in Sv (10^6/s)

         S_specific = [];
         theta_specific = [];
         V_specific = [];

	%compensation method below: might need to be changed for V2....
	for f = 1:cGrid %previously rGrid
	    %I am finding the discrepancy of the cumulative volume transfer at the bottom, then compensating by adding flow/number of cells with non-NaN values
	    %Q_cut = squeeze(Q(f,:,:));
        Q_cut = squeeze(Q(:,f,:));
	    Q_cut_copy = Q_cut; %copy for figuring out how many non NaN-s to divide among 
	    maskNaN = isnan(Q_cut_copy);
	    Q_cut_copy(~isnan(Q_cut_copy)) = 1;
	    Q_cut_copy(isnan(Q_cut_copy)) = 0;
	    sumReal = sum(sum(Q_cut_copy));

	    Q_cut_copy = Q_cut;
	    Q_cut_copy = nansum(Q_cut_copy,1);%changed from 2 to 1, considering changes in dimensions for input files
	    Q_cut_copy = cumsum(Q_cut_copy);
	    difference = Q_cut_copy(end); %Theoretically, this should be 0 but usually isn't due to resolution issues and average flows within a specific cell
	    offset = difference/sumReal;
	    Q_cut(~isnan(Q_cut)) = Q_cut(~isnan(Q_cut)) - offset; %offsets all non-NaN volumetric transfer cells along a certain latitude
	    %Q(f,:,:) = Q_cut; %corrects the original volumetric flow dataset by offset 
        Q(:,f,:) = Q_cut; %corrects the original volumetric flow dataset by offset 
	end

	Q_input_nc = cat(4,Q_input_nc,Q);
    %Q_cut = squeeze(Q(R,:,:)); %had a 3d array around 50N, now summing across just 50N
    Q_cut = squeeze(Q(:,R,:)); %had a 3d array around 50N, now summing across just 50N
    Q = [];

    %sigma0_cut = squeeze(sigma0_specific(R,:,:));
    sigma0_cut = squeeze(sigma0_specific(:,R,:));
	sigma0_input_nc = cat(4, sigma0_input_nc,sigma0_specific);

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
Q_filename = ['Q_V_gridded_' modelName '_instance' num2str(instance) '_' num2str(spacing) '_V2.nc'];
disp(Q_filename)
nccreate(string(Q_filename),"Flow", ...
"Dimensions",{"x",r_nc,"y", c_nc, "depth", z_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(Q_filename), "Flow", Q_input_nc)

Q_input_nc = [];
Q_filename = [];

sigma0_filename = ['Sigma0_gridded_' modelName '_instance' num2str(instance) '_' num2str(spacing) '_V2.nc'];
nccreate(string(sigma0_filename), "Sigma0", ...
"Dimensions", {"x", r_nc, "y", c_nc, "depth", z_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(sigma0_filename), "Sigma0", sigma0_input_nc)

sigma0_input_nc = [];
sigma0_filename = [];
%keeping .nc file extremely simple right now as I had previous problems with nccreate(), going to save other important details about the regridded file as a .mat file currently
end

%write .nc file for surface sigma0, sal, theta values
[r_nc,c_nc,t_nc] = size(sigma0_surface_nc);
sigma0_filename = ['Sigma0_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'];
nccreate(string(sigma0_filename), "Sigma0", "Dimensions", {"x", r_nc, "y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(sigma0_filename), "Sigma0", sigma0_surface_nc)

S_filename = ['S_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'];
nccreate(string(S_filename),"so", "Dimensions",{"x",r_nc,"y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(S_filename), "so", sal_surface_nc)

Theta_filename = ['Theta_surface_gridded_' modelName '_' num2str(spacing) '_V2.nc'];
nccreate(string(Theta_filename),"thetao", "Dimensions",{"x",r_nc,"y", c_nc, "time", t_nc}, "Format", "classic");
ncwrite(string(Theta_filename), "thetao", theta_surface_nc)

%% Calculating Annual Average MOC Strength and Average Sigma MOC, saving data
annual_max_Q_value = mean(reshape(max_Q_value, [12 length(max_Q_value)/12])); %assumes data is divisible by 12
meanSigmaMoc = mean(sigmaMocValue);

%make sure to rename these when starting a new model!
input1 = string(['sigma_and_MOC_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '_V2.mat']);
input2 = string(['regridded_dims_' modelName '_' num2str(spacing) '_V2.mat']);

save(input1, getVarName(time), getVarName(max_Q_value), getVarName(sigmaMocValue), ...
    getVarName(densitySpacing), getVarName(annual_max_Q_value)) 
save(input2, getVarName(time), getVarName(latitude), getVarName(longitude), ...
    getVarName(depth), getVarName(delta_z))

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
savefig(['SigmaMoc_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '_V2.fig'])

figure
plot(time, max_Q_value)
grid on
xlabel('Time [years]')
ylabel('Maximum Cumulative Flow (Sv) ')
title(['Maximum Cumulative Flow over 50N - Model: '  modelName ', Spacing: ' num2str(spacing)])
savefig(['max_Qflow_50N_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '_V2.fig'])

figure
plot(1850:2014, annual_max_Q_value, 'b', 'Linewidth', 1.5)
grid on
xlabel('Time [years]')
ylabel('Maximum Cumulative Flow (Sv) ')
title(['Annual Average Maximum Cumulative Flow over 50N - Model: '  modelName ', Spacing: ' num2str(spacing)])
savefig(['annual_avg_Qflow_50N_' modelName '_' num2str(spacing) '_ds_' num2str(densitySpacing) '_V2.fig'])


function out = getVarName(var)
    out = string(inputname(1));
end
