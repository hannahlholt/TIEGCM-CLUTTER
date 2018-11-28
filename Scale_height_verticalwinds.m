%   Author: Torfinn Johnsud
%   Created: Fall 2017  Modified: Aug 20 2018 by HLH
%   Plots scale heights for N2,O,He,and mean molecular mass. Also plots
%   vertical wind and total density profiles.
 

clear all
% close all
clc

aa1=['/home/haho3703/TIEGCM/TIEGCM_files/'];
linear=1;
alpha = -0.38;            % thermal diffusion factor

%----------------
ut_want=1;          % what time segment desired from simulation
feat = 4;           % Select latitude and longitude desired
pdrag = 1;          % 1 if pdrag file used, 0 if not

if feat == 1
lon_want=100;       % North nighttime feature
lat_want=63.75;     % North nighttime feature
picname='N_night';  % for graphing
end
if feat == 2
lon_want=102.5;     % South nighttime feature
lat_want=-58.75;    % South nighttime feature
picname='S_night';
end
if feat == 3
lat_want=58.75;     % Intermediate density area
lon_want=40;        % Intermediate density area
picname='Interm';
end
if feat == 4 
lon_want=-150;      % Daytime feature
lat_want=11.25;     % Daytime feature
picname='N_day';
end
if feat == 5
lon_want=-150;        % Select Longitude desired   
lat_want=-1.25;      % Select Latitude desired
picname='ETA_posCrestB';
end

%------------------
atom_unit=1.67e-27; % kg/unit
%-----------------

%----- Loading and Condensing Viki's TIEGCM simulation ---------------------
%---------------------------------------------------------------------------

% LOAD-------------------

% format is (lon,lat,ilev,UT) format
if pdrag == 1
    filename =[aa1, 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'];
    id = 'pdrag';
end
if pdrag == 0
    filename =[aa1, 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc'];
    id = 'no Ion Drag';
end

den = ncread(filename,'DEN');
zg = ncread(filename,'ZG');
he = ncread(filename,'HE');
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32);%Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');

% CONDENSE-----------------

% fixed UT time 
den_tp1=squeeze(den(:,:,:,ut_want+1));% Total neutral density
zg_tp1=squeeze(zg(:,:,:,ut_want+1));% Geometric height
% z_tp1=squeeze(z(:,:,:,ut_want+1));% Geopotential height
he_tp1=squeeze(he(:,:,:,ut_want+1));% Helium mass mixing ratio
N2_tp1=squeeze(n2(:,:,:,ut_want+1));
O1_tp1=squeeze(o1(:,:,:,ut_want+1));
tn_tp1=squeeze(tn(:,:,:,ut_want+1));% Neutral temperature
mbar_tp1=squeeze(mbar(:,:,:,ut_want+1));% Mean molecular mass
wn_tp1 = squeeze(wn(:,:,:,ut_want+1));


%fixed latitude and longitude  
i_index=find(lat==lat_want);

% Condense to selected latitude
den_alt_1=squeeze(den_tp1(:,i_index,:));
he_alt_1=squeeze(he_tp1(:,i_index,:));%TIEGCM output is mass mixing ratio
N2_alt_1=squeeze(N2_tp1(:,i_index,:));
O1_alt_1=squeeze(O1_tp1(:,i_index,:));
geometric_alt_1=squeeze(zg_tp1(:,i_index,:))/1e5; %Convert from cm to km
% geop_alt_1=squeeze(z_tp1(:,i_index,:))/1e5;
tn_alt_1=squeeze(tn_tp1(:,i_index,:));
mbar_alt_1=squeeze(mbar_tp1(:,i_index,:));
wn_alt_1 = squeeze(wn_tp1(:,i_index,:));

% -----Mean over the Lonitude-----
% den_alt=mean(den_alt_1);
% he_alt=mean(he_alt_1);
% N2_alt=mean(N2_alt_1);
% O1_alt=mean(O1_alt_1);
% % geop_alt=mean(geop_alt_1); 
% geom_alt=mean(geometric_alt_1);
% tn_alt=mean(tn_alt_1);
% mbar_alt=mean(mbar_alt_1);
% wn_alt = mean(wn_alt_1);

% -----Select the Longitude-----
i_index = find(lon==lon_want);
den_alt = squeeze(den_alt_1(i_index,:));
he_alt = squeeze(he_alt_1(i_index,:));
N2_alt = squeeze(N2_alt_1(i_index,:));
O1_alt = squeeze(O1_alt_1(i_index,:));
% geop_alt = squeeze(geop_alt_1(i_index,:));
geom_alt = squeeze(geometric_alt_1(i_index,:));
tn_alt = squeeze(tn_alt_1(i_index,:));
mbar_alt = squeeze(mbar_alt_1(i_index,:));%kg/kmol
wn_alt = squeeze(wn_alt_1(i_index,:));


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% Number Density
nhe_bf=he_alt.*den_alt*1e3/4/atom_unit;   %helium number density
nN2=N2_alt.*den_alt*1e3/28.01/atom_unit; % N2 number density
nO1=O1_alt.*den_alt*1e3/16/atom_unit; % O1 number density

%Mass Density
mhe = he_alt.*den_alt;% Helium Mass Density g/cm^3 
mN2 = N2_alt.*den_alt;% N2 Mass Density
mO1 = O1_alt.*den_alt;% Oxygen Mass Density


r=6372;%Earth Radius
g_tp=9.8.*r^2./(r+geom_alt).^2; %Find g as altitude changes
k=1.38e-23;% Boltzman's Constant
atom_unit=1.67e-27;%kg/unit
Av = 6.022141*10^23;
kmol=6.02214*10^26;
mmw_he=0.004; %Helium atomic mass
mmw_N2=0.02801; %N2 molecular mass
mmw_O1=0.016; % O1 molecular mass
points = length(tn_alt);
meanmass = mbar_alt/1000;


%-----Finding Pressure Scale Heights for gas constituents from kT/mg-----
Hp_he=k.*tn_alt./(mmw_he/Av.*g_tp)/1000;   %in Km
Hp_mean=k.*tn_alt./(meanmass/Av.*g_tp)/1000;
Hp_n2=k.*tn_alt./(mmw_N2/Av.*g_tp)/1000;
Hp_o1=k.*tn_alt./(mmw_O1/Av.*g_tp)/1000;

% -----Calculate Scale Heights Using 3 Point Differentiation Technique-----
H_he_star = zeros(points,1);
H_dentot = zeros(points,1);
H_n2 = zeros(points,1);
H_o1 = zeros(points,1);
H_tn = zeros(points,1);
H_mass = zeros(points,1);


for i=1:points          % ******* VALUES CALCULATED HERE ARE equal to 1/the  corrresponding scale heights!!!!!!
    if i==1 %First Point
        coeff1 = (2*geom_alt(1)-geom_alt(2)-geom_alt(3))/((geom_alt(1)-...
            geom_alt(2))*(geom_alt(1)-geom_alt(3)));
        
        coeff2 = (2*geom_alt(1)-geom_alt(1)-geom_alt(3))/((geom_alt(2)-...
            geom_alt(1))*(geom_alt(2)-geom_alt(3)));
        
        coeff3 = (2*geom_alt(1)-geom_alt(1)-geom_alt(2))/((geom_alt(3)-...
            geom_alt(1))*(geom_alt(3)-geom_alt(2)));
        
        H_he_star(1) = -1/mhe(1)*(mhe(1)*coeff1+mhe(2)*coeff2+...
            mhe(3)*coeff3);
        H_dentot(1) = -1/den_alt(1)*(den_alt(1)*coeff1+den_alt(2)*coeff2+...
            den_alt(3)*coeff3);
        H_n2(1) = -1/mN2(1)*(mN2(1)*coeff1+mN2(2)*coeff2+...
            mN2(3)*coeff3);
        H_o1(1) = -1/mO1(1)*(mO1(1)*coeff1+mO1(2)*coeff2+...
            mO1(3)*coeff3);
        H_tn(1) = 1/tn_alt(1)*(tn_alt(1)*coeff1+tn_alt(2)*coeff2+...
            tn_alt(3)*coeff3);
        H_mass(1) = -1/mbar_alt(1)*(mbar_alt(1)*coeff1+mbar_alt(2)*coeff2+...
            mbar_alt(3)*coeff3);
        
        
    elseif i==points %Last point
        coeff1 = (2*geom_alt(i)-geom_alt(i-1)-geom_alt(i))/((geom_alt(i-2)-...
            geom_alt(i-1))*(geom_alt(i-2)-geom_alt(i)));
        
        coeff2 = (2*geom_alt(i)-geom_alt(i-2)-geom_alt(i))/((geom_alt(i-1)-...
            geom_alt(i-2))*(geom_alt(i-1)-geom_alt(i)));
        
        coeff3 = (2*geom_alt(i)-geom_alt(i-2)-geom_alt(i-1))/((geom_alt(i)-...
            geom_alt(i-2))*(geom_alt(i)-geom_alt(i-1)));
        
        H_he_star(i) = -1/mhe(i)*(mhe(i-2)*coeff1+mhe(i-1)*coeff2+...
            mhe(i)*coeff3);
        H_dentot(i) = -1/den_alt(i)*(den_alt(i-2)*coeff1+den_alt(i-1)*coeff2+...
            den_alt(i)*coeff3);
        H_n2(i) = -1/mN2(i)*(mN2(i-2)*coeff1+mN2(i-1)*coeff2+...
            mN2(i)*coeff3);
        H_o1(i) = -1/mO1(i)*(mO1(i-2)*coeff1+mO1(i-1)*coeff2+...
            mO1(i)*coeff3);
        H_tn(i) = 1/tn_alt(i)*(tn_alt(i-2)*coeff1+tn_alt(i-1)*coeff2+...
            tn_alt(i)*coeff3);
        H_mass(i) = -1/mbar_alt(i)*(mbar_alt(i-2)*coeff1+mbar_alt(i-1)*coeff2+...
            mbar_alt(i)*coeff3);

        
    else %Middle Points
        coeff1 = (2*geom_alt(i)-geom_alt(i)-geom_alt(i+1))/((geom_alt(i-1)-...
            geom_alt(i))*(geom_alt(i-1)-geom_alt(i+1)));
        
        coeff2 = (2*geom_alt(i)-geom_alt(i-1)-geom_alt(i+1))/((geom_alt(i)-...
            geom_alt(i-1))*(geom_alt(i)-geom_alt(i+1)));
        
        coeff3 = (2*geom_alt(i)-geom_alt(i-1)-geom_alt(i))/((geom_alt(i+1)-...
            geom_alt(i-1))*(geom_alt(i+1)-geom_alt(i)));
        
        H_he_star(i) = -1/mhe(i)*(mhe(i-1)*coeff1+mhe(i)*coeff2+...
            mhe(i+1)*coeff3);
        H_dentot(i) = -1/den_alt(i)*(den_alt(i-1)*coeff1+den_alt(i)*coeff2+...
            den_alt(i+1)*coeff3);
        H_n2(i) = -1/mN2(i)*(mN2(i-1)*coeff1+mN2(i)*coeff2+...
            mN2(i+1)*coeff3);
        H_o1(i) = -1/mO1(i)*(mO1(i-1)*coeff1+mO1(i)*coeff2+...
            mO1(i+1)*coeff3);
        H_tn(i) = 1/tn_alt(i)*(tn_alt(i-1)*coeff1+tn_alt(i)*coeff2+...
            tn_alt(i+1)*coeff3);
        H_mass(i) = -1/mbar_alt(i)*(mbar_alt(i-1)*coeff1+mbar_alt(i)*coeff2+...
            mbar_alt(i+1)*coeff3);

    end
end

% ****** Now the ACTUAL scale heights are found below....

%-----Get Scale Height from Inverse-----
H_he_star = 1./H_he_star;
H_den_tot = 1./H_dentot;
H_N2_star = 1./H_n2;
H_O1_star = 1./H_o1(4:end);
H_temp = (1./H_tn).';
H_temp_he = H_temp/(1 + alpha);             % thermal diffusion factor alpha = -0.38
H_mass = (1./H_mass).';

if linear==1
    %--------Linear Gradient Calculations--------
    % -----Helium Scale Height-----
    lnn_bf = reallog(nhe_bf);
    H_he_star_r_bf = zeros(points,1);
    for i=1:points
        if i==1
            H_he_star_r_bf(i)=-(lnn_bf(2)-lnn_bf(1))./...
                            (geom_alt(2)-geom_alt(1));
        elseif i==points
            H_he_star_r_bf(i)=-(lnn_bf(i)-lnn_bf(i-1))./...
                            (geom_alt(i)-geom_alt(i-1));
        else
            H_he_star_r_bf(i)=-(lnn_bf(i+1)-lnn_bf(i-1))./...
                            (geom_alt(i+1)-geom_alt(i-1));
        end
    end

    H_he_star=1./H_he_star_r_bf;

    % -----Scale Height for Total Density-----
    lnnden = reallog(den_alt);
    H_dentot = zeros(points,1);
    for i=1:points
        if i==1
            H_dentot(1) = -(lnnden(2)-lnnden(1))/(geom_alt(2)-...
                geom_alt(1));
        elseif i==points
            H_dentot(i) = -(lnnden(i)-lnnden(i-1))/(geom_alt(i)-...
                geom_alt(i-1));
        else
        H_dentot(i) = -(.5*(lnnden(i)-lnnden(i-1))/(geom_alt(i)-...
            geom_alt(i-1))+.5*(lnnden(i+1)-lnnden(i))/(geom_alt(i+1)-...
            geom_alt(i)));
        end
    end
    H_den_tot = 1./H_dentot;

    % -----Scale Height for N2-----
    lnnN2 = reallog(mN2);
    H_n2 = zeros(points,1);
    for i=1:points
        if i==1
            H_n2(1) = -(lnnN2(2)-lnnN2(1))/(geom_alt(2)-...
                geom_alt(1));
        elseif i==points
            H_n2(i) = -(lnnN2(i)-lnnN2(i-1))/(geom_alt(i)-...
                geom_alt(i-1));
        else
        H_n2(i) = -(lnnN2(i+1)-lnnN2(i-1))/(geom_alt(i+1)-...
            geom_alt(i-1));
        end
    end
    H_N2_star = 1./H_n2;

    % -----Scale Height for O1-----
    lnO1 = reallog(nO1);
    H_o1 = zeros(points,1);
    for i=1:points
        if i==1
            H_o1(1) = -(lnO1(2)-lnO1(1))/(geom_alt(2)-...
                geom_alt(1));
        elseif i==points
            H_o1(i) = -(lnO1(i)-lnO1(i-1))/(geom_alt(i)-...
                geom_alt(i-1));
        else
        H_o1(i) = -(lnO1(i+1)-lnO1(i-1))/(geom_alt(i+1)-...
            geom_alt(i-1));
        end
    end
    H_O1_star = 1./H_o1(4:end);
end


%----Put Together Diffusive Profiles and Mean Mass Profile-----
H_he_diff = 1./(1./H_temp_he+1./Hp_he);% Helium diffusive profile
H_he_diff_NO_thermal = 1./(1./H_temp+1./Hp_he); % helium diffusive profile with no thermal difffusion
H_tot = 1./(1./H_temp+1./Hp_mean+1./H_mass);% Mean mass profile
H_N2_diff = 1./(1./H_temp+1./Hp_n2);% N2 Diffusive profile
H_O1_diff = 1./(1./H_temp+1./Hp_o1);% O1 Diffusive profile

%H_he_difference_T = (H_he_diff - H_he_diff_NO_thermal) ./ H_he_diff .* 100;        % percent difference of helium diffusive profiles with and w/o thermal diffusion
H_he_difference_tot = zeros(length(H_he_star),1);

for i=1:length(H_he_star)
   H_he_difference_tot(i) = (H_he_star(i) - H_he_diff(i)) / H_he_star(i) * 100; % percent difference between TIEGCM helium scale height and pure diffusive profile (w/ thermal diff.)
   
end
                 
% -----Plot the Scale Heights-----

%--- HELIUM
start_pos = 10;
end_pos = 47;

f=figure;
set(f, 'position', [0, 500, 650, 900])
subplot(2,1,1)
hold on;
plot(H_he_star(start_pos:end_pos), geom_alt(start_pos: end_pos),'r','linewidth', 2);
%plot(H_he_star_ln(1:end-2), geom_alt(1:end-2));
plot(H_he_diff(start_pos:end_pos),geom_alt(start_pos:end_pos),'b','linewidth', 2);
plot(H_tot(start_pos:end_pos),geom_alt(start_pos:end_pos),'color',[.2,.7,.2],'linewidth', 2);
legend('H_i^* TIEGCM','H_i Diffusive Eq','H Total' ,'location','north');
title(['Helium Scale Heights, lat/lon=',num2str(lat_want),'/',num2str(lon_want),', UT=0.03, ',id, ', alpha =', num2str(alpha)]);
xlabel('Scale Height (km)');
ylabel('Geometric Altitude (km)');
grid on;
% % saveas(f,'H_he_intermediate_density.png');
% 
subplot(2,1,2) 
plot(H_he_difference_tot(start_pos:end_pos), geom_alt(start_pos:end_pos),'b','linewidth', 2);
%legend('Therm Diffusion - No Therm Diffusion')
title('TIEGCM - Diffusive Equilibrium Difference');
xlabel('Percent Difference (%)');
ylabel('Geometric Altitude (km)');
name = ['~/Desktop/He_diffusive_difference_perc_', picname, '_pdrag.png'];
grid on;
% saveas(f, name);

%---- ATOMIC OXYGEN
% f=figure;
% hold on
% plot(H_O1_star(3:44),geom_alt(6:47),'r','linewidth', 2);
% % plot(H_o1_star_ln(6:end-2), geom_alt(6:end-2));
% plot(H_O1_diff(1:47),geom_alt(1:47),'b','linewidth', 2);
% plot(H_tot(1:47),geom_alt(1:47),'color',[.2,.7,.2],'linewidth', 2);
% legend('H_i^* TIEGCM','H_i Diffusive Eq','H Total' ,'location','north');
% title(['Oxygen Scale Heights lon=',num2str(lon_want),' lat=',num2str(lat_want),' UT=0.03 ',id]);
% xlabel('Scale Height (km)');
% ylabel('Geometric Altitude (km)');
% % saveas(f,'H_O1_intermediate_density.png');
% 

%---- MOLECULAR NITROGEN
% f=figure;
% hold on
% plot(H_N2_star(1:47),geom_alt(1:47),'r','linewidth', 2);
% % plot(H_n2_star_ln(1:end-2), geom_alt(1:end-2));
% plot(H_N2_diff(1:47),geom_alt(1:47),'b','linewidth', 2);
% plot(H_tot(1:47),geom_alt(1:47),'color',[.2,.7,.2],'linewidth', 2);
% legend('H_i^* TIEGCM','H_i Diffusive Eq','H Total' ,'location','northwest');
% title(['N2 Scale Heights lon=',num2str(lon_want),' lat=',num2str(lat_want),' UT=0.03 ',id]);
% xlabel('Scale Height (km)');
% ylabel('Geometric Altitude (km)');
% % saveas(f,'H_N2_intermediate_density.png');
% 

%----- COMBO
% figure
% hold on
% plot(H_he_star(1:47), geom_alt(1:47),'r','linewidth', 2);
% plot(H_O1_star(3:44),geom_alt(6:47),'b','linewidth', 2);
% plot(H_N2_star(1:47),geom_alt(1:47),'color',[.2,.7,.2],'linewidth', 2);
% legend('H_H_e^*','H_O_1^*','H_N_2^*','location','north');
% xlabel('Scale Height (km)');
% ylabel('Geometric Altitude (km)');
% title(['Scale Heights lon=',num2str(lon_want),' lat=',num2str(lat_want),' UT=0.03 ',id])

% figure
% plot (reallog(mhe),geom_alt);

% ------Mean Mass of Thermosphere
% plot(meanmass/Av/atom_unit, geom_alt, 'linewidth', 2);
% title(['Mean Mass of Thermosphere, lat/lon=',num2str(lat_want),'/',num2str(lon_want),', UT=0.03, ',id]);
% xlabel('Mean mass [amu]');
% ylabel('Altitude [km]');
% grid on;

% 
% figure
% subplot(2,1,1)
% hold on;
% plot(H_temp(1:30), geom_alt(1:30),'linewidth',2);
% plot(Hp_mean(1:end-1), geom_alt(1:end-1),'linewidth',2);
% plot(H_mass(:),geom_alt(:),'linewidth',2);
% plot(H_tot(1:50),geom_alt(1:50),'color',[.2,.7,.2],'linewidth', 2);
% legend('Temp','Pressure','Mean mass', 'Total');
% xlabel('Scale Height [km]');
% ylabel('Geometric Altitude [km]');
% title(['Scale Heights lat/lon=',num2str(lat_want),'/',num2str(lon_want),' UT=0.03, ',id]);
% grid on;
% 
% subplot(2,1,2)
% plot(Hp_mean(1:end-1), geom_alt(1:end-1),'color','r','linewidth', 2);
% legend('Pressure', 'location', 'northwest');
% xlabel('Scale Height [km]');
% ylabel('Geometric Altitude [km]');
% grid on;

% 
% figure
% plot(meanmass,geom_alt);
% ylabel('Geometric Altitude [km]');
% xlabel('Mean Molecular Mass [kg/kmol]');
% title(['Mean Molecular Mass vs Altitude lon=',num2str(lon_want),' lat=',num2str(lat_want),' UT=0.03 no Ion Drag'])
% 
% %-----Plot Density and Vertical winds for mass flux-----
% figure
% plot(wn_alt,geom_alt,'b','linewidth',2);
% xlabel('Vertical Winds (cm/s)');
% ylabel('Geometric Altitude (km)');
% title(['Vertical Winds With Altitude, lat/lon=',num2str(lat_want),'/',num2str(lon_want),' UT=0.03, ',id]);
% grid on;
% 
% figure
% plot(den_alt,geom_alt,'r','linewidth',2);
% set(gca, 'Xscale', 'log');
% xlabel('Log Total Density (g/cm^3)');
% ylabel('Geometric Altitude (km)');
% title('Total Mass Density');

% wn_alt_day = wn_alt;
% wn_alt_night_south = wn_alt;
% wn_alt_night_north = wn_alt;
% wn_alt_intermediate = wn_alt;

% den_alt_day = den_alt;
% den_alt_night_south = den_alt;
% den_alt_night_north = den_alt;
% den_alt_intermediate = den_alt;

% figure
% plot(H_temp(1:end-30),geom_alt(1:end-30));
% xlabel('Temperature Scale Height [km]');
% ylabel('Altitude [km]');
% title('Nighttime Feature');

% figure
% plot(Hp_o1(1:end-30),geom_alt(1:end-30));
% xlabel('Pressure Scale Height [km]');
% ylabel('Altitude [km]');
