function [  ] = Massfluxplots(wn_day,wn_nn,wn_ns,wn_i,den_day,den_nn,den_ns,den_i,altitude)
%Plots vertical winds and density with altitude. Mass flux is
%density*wind.

figure
hold on
plot(wn_day,altitude,'r','linewidth',2);
plot(wn_i,altitude,'b','linewidth',2);
plot(wn_nn,altitude,'k','linewidth',2);
plot(wn_ns,altitude,'--k','linewidth',2);
legend('Day Feature','Intermediate Location','Night North Feature','Night South Feature','location','southwest');
xlabel('Vertical Wind (cm/s)');
ylabel('Geometric Altitude (km)');
title('Vertical Winds with Altitude');

figure 
hold on
plot(den_day,altitude,'r');
plot(den_i,altitude,'b');
plot(den_nn,altitude,'k');
plot(den_ns,altitude,'--k');
set(gca,'Xscale','log');
legend('Day Feature','Intermediate Location','Night North Feature','Night South Feature','location','southwest');
xlabel('Mass Density (g/cm^3)');
ylabel('Geometric Altitude (km)');
title(' Total Mass Density with Altitude');
end

