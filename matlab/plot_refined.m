tau = [0.01, 0.08, 0.4, 1.0];

unpinned = [0.501     ,0.570    ,0.890     ,1.479];
normal_4 = [0.48188   ,0.555    ,0.874     ,1.4659];
normal_8 = [0.4532    ,0.525    ,0.842     ,1.451];
normal_12 = [0.4503   ,0.519    ,0.8364    ,1.42913];
normal_16 = [0.431493 ,0.4994   ,0.814747  ,1.41275];
normal_20 = [0.410498 ,0.480    ,0.796215  ,1.39092];
normal_24 = [0.391599 ,0.461333 ,0.774648  ,1.37455];
normal_28 = [0.36664  ,0.433778 ,0.750836  ,1.34782];
normal_32 = [0.349648 ,0.419573 ,0.741193  ,1.34209];

sig_1 = [0.501, 0.48188, 0.4532, 0.4503, 0.431493, 0.410498, 0.391599, 0.36664, 0.349648] / 0.501;
sig_2 = [0.570, 0.555, 0.525, 0.519, 0.4994, 0.480, 0.461333, 0.433778, 0.419573] / 0.570;
sig_3 = [0.890, 0.874, 0.842, 0.8364, 0.814747, 0.796215, 0.774648, 0.750836, 0.741193] / 0.890;
sig_4 = [1.479,1.4659,1.451,1.42913,1.41275,1.39092,1.37455,1.34782,1.34209] / 1.479;
percentages = [0.0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32];

figure(1)
plot(percentages, sig_1,'Marker','o')
hold on
plot(percentages, sig_2,'Marker','o')
hold on
plot(percentages, sig_3,'Marker','o')
hold on
plot(percentages, sig_4,'Marker','o')
hold on
grid on
xlabel(" Fraction  of pinning")
ylabel("\sigma/\sigma_{0}")
grid on
legend('\tau=0.01 ','\tau=0.08','\tau=0.4','\tau=1.0');

figure(2)
semilogx(tau,unpinned,'Marker','s');
hold on
semilogx(tau,normal_4,'Marker','s');
hold on
semilogx(tau,normal_8,'Marker','s');
hold on
semilogx(tau,normal_12,'Marker','s');
hold on
semilogx(tau,normal_16,'Marker','s');
hold on
semilogx(tau,normal_20,'Marker','s');
hold on
semilogx(tau,normal_24,'Marker','s');
hold on
semilogx(tau,normal_28,'Marker','s');
hold on
semilogx(tau,normal_32,'Marker','s');
xlabel("\tau (k_BT/a_0^2)")
ylabel("\sigma (k_BT/a_0^2)")
grid on
legend('unpinned','4%','8%','12%','16%','20%','24%','28%','32%');
