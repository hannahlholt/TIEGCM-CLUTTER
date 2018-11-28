%This was used to produce the results in the error analysis of the 
%exponential gradient technique. Data from 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc'
%and the scale height code was used. 

format long
i = 15;
x = linspace(geom_alt(i-2),geom_alt(i+2),200);
plot(x,mhe(i)*exp((x-geom_alt(i))/-H_he_star(i)),'color', 'k')
hold on
plot(geom_alt, mhe, 'color', 'r');
xlabel('Altitude (km)');
ylabel('Density');
title('Approximation Curve and Data');
scatter([geom_alt(i-1), geom_alt(i+1)],[mhe(i-1), mhe(i+1)],'MarkerEdgeColor','b');
legend('Approximation','Data','Adjacent Data Points');

% z = mhe(36)/(mhe(35)*exp((geom_alt(36)-geom_alt(35))/-H_he_star_ln(35)));
% y = mhe(34)/(mhe(35)*exp((geom_alt(34)-geom_alt(35))/-H_he_star_ln(35)));
% % err1 = log(z/y)/(geom_alt(36)-geom_alt(34))
% err2 = (H_he_star(35)*(geom_alt(36)-geom_alt(35))-H_he_star(35)*(geom_alt(34)-geom_alt(35))+log(mhe(36)/mhe(34)))/(geom_alt(36)-geom_alt(34))
% N = (geom_alt(36)-geom_alt(35))/-H_he_star_ln(35);
% M = (geom_alt(34)-geom_alt(35))/-H_he_star_ln(35);
% err3 = (1/M-1/N+log(mhe(36)/mhe(34)))/(geom_alt(36)-geom_alt(34))

% top = log(mhe(i-1))-log(mhe(i+1))+(geom_alt(i+1)-geom_alt(i-1))/H_he_star(i);
% err = top/(geom_alt(i+1)-geom_alt(i-1));
lambda = 1/H_he_star(i);
% percentage = err/lambda

y1 = mhe(i-1)*10^15;
y3 = mhe(i+1)*10^15;
fx1 = mhe(i)*exp((geom_alt(i-1)-geom_alt(i))/-H_he_star(i))*10^15;
fx3 = mhe(i)*exp((geom_alt(i+1)-geom_alt(i))/-H_he_star(i))*10^15;

eta = y1/fx1;
delta = y3/fx3;

err = log(eta/delta)/(geom_alt(i-1)-geom_alt(i+1));
percentage = err/lambda*100;

err1 = H_he_star(i)*log(eta/delta)