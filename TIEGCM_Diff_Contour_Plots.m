%Author: Torfinn Johnsrud
%Created: Fall 2017  UPDATED FOR SPECIFIC ANALYSIS 8/23/2018 by HLH - see
%Torfinn backop code for original file
% This program blends the Scale_height verticalwinds.m with the orginal
% TIEGCM_contour.m program. 

clear all;
close all;
clc;

aa1 = '/home/haho3703/TIEGCM/TIEGCM_files/';
aa2 = '/home/haho3703/TIEGCM/Contour_textfiles/';
linear=1;

%----------------
ut_want = 1;%
alt_want = 400;
z0 = 120;                   % bottom integration value [km]
pdrag = 1;


%-----Loading Viki's tiegcm simulation-----
% Follows (lon,lat,ilev,UT) format
if pdrag == 1
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename = [aa1, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'ctrSS';
end

den = ncread(filename,'DEN');           % [g/cm^3]
zg = ncread(filename,'ZG')/1e5;         % geometric height [km]
he = ncread(filename,'HE');
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32);%Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');

% fixed UT time 
den=squeeze(den(:,:,:,ut_want+1));      % Total neutral density [g/cm^3]
zg_tp1=squeeze(zg(:,:,:,ut_want+1));    % Geometric height [km]
% z_tp1=squeeze(z(:,:,:,ut_want+1));    % Geopotential height
he=squeeze(he(:,:,:,ut_want+1));        % Helium mass mixing ratio
tn=squeeze(tn(:,:,:,ut_want+1));        % Neutral temperature [K]
mbar=squeeze(mbar(:,:,:,ut_want+1));    % Mean molecular mass [kg/kmol]

% now each matrix above has a value for every lat, lon, and pressure level. We
% care about helium, temperature, and mean moleulcar mass. 
%-------------------------------------------------------------------------

r=6372;                                 % Earth Radius [km]
g_tp=9.8.*r^2./(r+zg_tp1).^2;           % Find g as altitude changes
k=1.38e-23;                             % Boltzman's Constant
atom_unit=1.67e-27;                     % [kg/unit]
Av = 6.022141*10^23;                    % [#/mol]
kmol=6.02214*10^26;
mmw_he=0.004;                           % Helium atomic mass [kg/mol]
mmw_N2=0.02801;                         % N2 molecular mass
mmw_O1=0.016;                           % O1 molecular mass
meanmass = mbar./1000;                  % mean molecular mass [kg/mol]

% ---- Helium Number Density
nhe_bf = (he.*den).*Av./(mmw_he.*1000); % actual helium number density [#/cm^3]

%------ Helium Mass Density
mhe = he.*den;                          % Helium Mass Density [g/cm^3] 

%-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
Hp_he=k.*tn./(mmw_he/Av.*g_tp)/1000;   %in Km
Hp_mean=k.*tn./(meanmass/Av.*g_tp)/1000;


%----- Finding temperature scale height for Helium -----------------------
%-------------------------------------------------------------------------
% Calculate using 3 Point Differentiation Technique

% 144 longitudes (360/2.5), 72 latitudes (180/2.5), 57 altitudes
lon_num = 144;
lat_num = 72;
alt_num = 57;
H_he_star = zeros(size(mhe));
H_tn = zeros(size(tn));
nhe_diff = zeros(lon_num, lat_num, alt_num);        % number density of helium at z if pure diff. profile.

% linearly interpolate to get these ...
nhe_real_120 = zeros(lon_num, lat_num, 1);          % actual number density of helium at z = 120 km
nhe_real_400 = zeros(lon_num, lat_num, 1);          % actual number density of helium at z = 400 km
nhe_diff_400 = zeros(lon_num, lat_num, 1);          % calculated number density of helium at z = 400 km

           
% ******* VALUES CALCULATED HERE ARE equal to 1/the  corrresponding scale heights!!!!!!

for i = 1:lon_num
    for j = 1:lat_num
        for k = 1:alt_num          
            if k == 1 % First Point   
                coeff1 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 2)-zg_tp1(i, j, 3))/((zg_tp1(i, j, 1)-...
                    zg_tp1(i, j, 2))*(zg_tp1(i, j, 1)-zg_tp1(i, j, 3)));

                coeff2 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 1)-zg_tp1(i, j, 3))/((zg_tp1(i, j, 2)-...
                    zg_tp1(i, j, 1))*(zg_tp1(i, j, 2)-zg_tp1(i, j, 3)));

                coeff3 = (2*zg_tp1(i, j, 1)-zg_tp1(i, j, 1)-zg_tp1(i, j, 2))/((zg_tp1(i, j, 3)-...
                    zg_tp1(i, j, 1))*(zg_tp1(i, j, 3)-zg_tp1(i, j, 2)));

                H_he_star(i,j, 1) = -1/mhe(i, j, 1)*(mhe(i, j, 1)*coeff1+mhe(i, j, 2)*coeff2+...
                    mhe(i, j, 3)*coeff3);
                H_tn(i, j, 1) = 1/tn(i, j, 1)*(tn(i, j, 1)*coeff1+tn(i, j, 2)*coeff2+...
                    tn(i, j, 3)*coeff3);


            elseif k == alt_num %Last point
                coeff1 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k))/((zg_tp1(i, j, k-2)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k-2)-zg_tp1(i, j, k)));

                coeff2 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-2)-zg_tp1(i, j, k))/((zg_tp1(i, j, k-1)-...
                    zg_tp1(i, j, k-2))*(zg_tp1(i, j, k-1)-zg_tp1(i, j, k)));

                coeff3 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-2)-zg_tp1(i, j, k-1))/((zg_tp1(i, j, k)-...
                    zg_tp1(i, j, k-2))*(zg_tp1(i, j, k)-zg_tp1(i, j, k-1)));

                H_he_star(i, j, k) = -1/mhe(i, j, k)*(mhe(i, j, k-2)*coeff1+mhe(i, j, k-1)*coeff2+...
                    mhe(i, j, k)*coeff3);
                H_tn(i, j, k) = 1/tn(i, j, k)*(tn(i, j, k-2)*coeff1+tn(i, j, k-1)*coeff2+...
                    tn(i, j, k)*coeff3);


            else %Middle Points
                coeff1 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k)-zg_tp1(i, j, k+1))/((zg_tp1(i, j, k-1)-...
                    zg_tp1(i, j, k))*(zg_tp1(i, j, k-1)-zg_tp1(i, j, k+1)));

                coeff2 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k+1))/((zg_tp1(i, j, k)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k)-zg_tp1(i, j, k+1)));

                coeff3 = (2*zg_tp1(i, j, k)-zg_tp1(i, j, k-1)-zg_tp1(i, j, k))/((zg_tp1(i, j, k+1)-...
                    zg_tp1(i, j, k-1))*(zg_tp1(i, j, k+1)-zg_tp1(i, j, k)));

                H_he_star(i, j, k) = -1/mhe(i, j, k)*(mhe(i, j, k-1)*coeff1+mhe(i, j, k)*coeff2+...
                    mhe(i, j, k+1)*coeff3);
                H_tn(i, j, k) = 1/tn(i, j, k)*(tn(i, j, k-1)*coeff1+tn(i, j, k)*coeff2+...
                    tn(i, j, k+1)*coeff3);

            end
        end
    end
end

% ****** Now the ACTUAL scale heights are found....

%-----Get Scale Height from Inverse-----
% The values below are the scale heights for every lat, long and altitude 

alpha = -0.38;                           % thermal diffusion factor
H_he_star = 1./H_he_star;                % actual scale height from density fit
H_temp = (1./H_tn);
H_temp_he = H_temp./(1 + alpha);         % helium temperature scale height
H_he_diff = 1./(1./H_temp_he+1./Hp_he);  % helium diffusive profile (i.e. the diffusive density scale height)

% -----Select Altitude-----   
geom_index_1 = zg_tp1-z0;           % starting altitude
geom_index_2 = zg_tp1-alt_want;     % ending altitude

for m = 1:lon_num 
    for n = 1:lat_num
            geom_index_start=geom_index_1(m,n,:);
            geom_index_start=squeeze(geom_index_start);
            [val_1, i_index_1] = min(abs(geom_index_start));    % Find altitude index closest to starting point
            
            geom_index_end=geom_index_2(m,n,:);
            geom_index_end=squeeze(geom_index_end);  
            [val_2, i_index_2] = min(abs(geom_index_end));      % Find altitude index to closest end point 
        
            for h_prime = i_index_1:(i_index_2+1)
                sum = 0;
                for h = i_index_1:h_prime
                    dz = zg_tp1(m, n, h+1) - zg_tp1(m, n, h);                       % altitude stepsize
                    avg = mean([1/H_he_diff(m, n, h+1), 1/H_he_diff(m, n, h)]);     % find average between scale height points
                    sum = sum + (avg * dz);                                         % integrate 
                end
                
                % this has non-zero values for every altitude between z0
                % and alt_want
                nhe_diff(m, n, h_prime) = nhe_bf(m, n, i_index_1) * exp(-sum);    % calculate density from exp. of scale height integral    
            end
            
            
            % -------------------------------------------------------------
            % -------------------------------------------------------------
            % Since the TIEGCM output is on pressure levels and not
            % altitudes, we can to interpolate the densities in between the
            % values to get the altitudes we want. We have to do this for
            % both the start and the end of the real He densities   
      
            
            % STARTING VALUES at Z ~ 120 km
            if zg_tp1(m,n,i_index_1) >= z0 % Linearly interpolate backwards
                x0 = zg_tp1(m,n,i_index_1-1);
                x1 = zg_tp1(m,n,i_index_1);
                y0 = nhe_bf(m,n,i_index_1-1);
                y1 = nhe_bf(m,n,i_index_1);
                
                nhe_real_120(m, n, 1) = (y0*(x1-z0)+y1*(z0-x0))/(x1-x0);
            end
            
            if zg_tp1(m,n,i_index_1) < z0 % Linearly interpolate forwards
                x0 = zg_tp1(m,n,i_index_1);
                x1 = zg_tp1(m,n,i_index_1+1);
                y0 = nhe_bf(m,n,i_index_1);
                y1 = nhe_bf(m,n,i_index_1+1);
                
                nhe_real_120(m, n, 1) = (y0*(x1-z0)+y1*(z0-x0))/(x1-x0);
            end

            % ENDING VALUES at Z ~ 400 km
            if zg_tp1(m,n,i_index_2) >= alt_want % Linearly interpolate backwards
                x0 = zg_tp1(m,n,i_index_2-1);
                x1 = zg_tp1(m,n,i_index_2);
                
                y0 = nhe_diff(m,n,i_index_2-1);
                y1 = nhe_diff(m,n,i_index_2);  
                nhe_diff_400(m, n, 1) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);             
                
                y0 = nhe_bf(m,n,i_index_2-1);
                y1 = nhe_bf(m,n,i_index_2);  
                nhe_real_400(m, n, 1) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);                
            end
            
            if zg_tp1(m,n,i_index_2) < alt_want % Linearly interpolate forwards
                x0 = zg_tp1(m,n,i_index_2);
                x1 = zg_tp1(m,n,i_index_2+1);
                
                y0 = nhe_diff(m,n,i_index_2);
                y1 = nhe_diff(m,n,i_index_2+1);
                nhe_diff_400(m, n, 1) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
                
                y0 = nhe_bf(m,n,i_index_2);
                y1 = nhe_bf(m,n,i_index_2+1);
                nhe_real_400(m, n, 1) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
            end
            % -------------------------------------------------------------
            % -------------------------------------------------------------                 
            
    end
end

% now we want to compare the actual density at ~400km to the calcualted
% density at ~400km
%------------------------------------------------------------------------

% find percent difference of diffusive - actual densities @ every lat/lon
perc_diff = ((nhe_real_400 ./ nhe_diff_400) - 1) .* 100;
perc_diff_OUTPUT = perc_diff.';                  % turns lon/lat ---> lat/lon
dlmwrite([aa2, 'He_dens_RealVsDiff_', id, '.txt'], perc_diff_OUTPUT);

% to look at the calculated diffusive He density
nhe_diff_OUTPUT = nhe_diff_400.';
dlmwrite([aa2, 'He_dens_Calculated_Diff_', id, '.txt'], nhe_diff_OUTPUT);

% to look at the real model helium density
nhe_real_OUTPUT = nhe_real_400.';
dlmwrite([aa2, 'He_dens_model_400km_', id, '.txt'], nhe_real_OUTPUT);

% to look at the real model helium density
nhe_real_120_OUTPUT = nhe_real_120.';
dlmwrite([aa2, 'He_dens_model_120km_', id, '.txt'], nhe_real_120_OUTPUT);





