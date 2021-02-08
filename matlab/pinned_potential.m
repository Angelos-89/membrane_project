std = 0.1; % mesured in a_0
hmin = -std/2; 
hmax = std/2;
dh = 1e-3;
h = hmin:dh:hmax;

kk = 1000;
E = 0.5*kk*h.*h; % E measured in k_bT

figure(1)
plot(h,E);
title("Potential energy of pinned site");
xlabel("h (a_0)");
ylabel("Energy (k_bT)");
ylim([0,1.5])
parallel = ones(length(h));
line(h, parallel, 'Color', 'Red');
grid on