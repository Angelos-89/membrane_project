std = 0.75512646208; % value obtained by a simulation with frame tension tau = 0
pinned_std = std/64; 

hmin = -1.5*pinned_std; 
hmax = +1.5*pinned_std;
dh = 1e-3;
h = hmin:dh:hmax;

kk = 2/(pinned_std^2); % demand the potential energy to be equal to k_bT = 1 (in dimensionless units) [k] = k_bT/a_0^2
E = 0.5*kk*h.*h;

figure(1)
plot(h,E);
title("Potential energy of pinned site");
xlabel("Height (a_0)");
ylabel("Energy (k_bT)");
xlim([hmin,hmax])
ylim([0,1.5])
parallel = ones(length(h));
vertical = ones(length(h));

yline(1,'r',"k_bT=1")
xline(-pinned_std,"k--","-std");
xline(+pinned_std,"k--","+std");
text(-0.003,0.3,"k \approx 14366")
grid on