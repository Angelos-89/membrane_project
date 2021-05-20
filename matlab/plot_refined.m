%----------------------------sigma-tau plot-----------------------%
tau_shiba = [0.01, 0.02, 0.04, ...
    0.1, 0.2, 0.4, ...
    1.0, 2.0, 4.0, ...
    10, 40, ...
    100, ...
    1000];
   
sigma_shiba = [0.509713537, 0.519494765, 0.538834329, ...
    0.597235415, 0.694884157, 0.890458633, ...
    1.47868855, 2.462251327, 4.435035741, ...
    10.37539235, 40.23293487, ...
    100.1352802, ...
    1000.018997]; 

tau = [0.01, 0.02, 0.04, 0.08, ...
       0.1, 0.2, 0.4, 0.8, ...
       1.0, 2.0, 4.0, 8.0, ...
       10.0,40.0...
       100.0, ...
       1000.0];
   
sigma_shiba = [0.509713537, 0.519494765, 0.538834329, ...
    0.597235415, 0.694884157, 0.890458633, ...
    1.47868855, 2.462251327, 4.435035741, ...
    10.37539235, 40.23293487, ...
    100.1352802, ...
    1000.018997]; 

tau = [0.01, 0.02, 0.04, 0.08, ...
       0.1, 0.2, 0.4, 0.8, ...
       1.0, 2.0, 4.0, 8.0, ...
       10.0,40.0...
       100.0, ...
       1000.0];
   
sigma_unpinned = [0.506143, 0.515, 0.5348571, 0.5728571 ...
         0.5938709, 0.69, 0.89, 1.2805 ...
         1.479, 2.465, 4.448, 8.41 ...
         10.4, 40.280, 100.18, 1000.03];
            
sigma_normal_8 = [0.46928, 0.47666, 0.4972, 0.5359677 ...
         0.5594, 0.653, 0.85666, 1.2496825 ...
         1.451, 2.420, 4.42, 8.39066 ...
         10.379677, 40.264, 100.17, 1000.03];

sigma_normal_16 = [0.427, 0.438, 0.4602, 0.496 ...
         0.519333, 0.615, 0.816, 1.21333 ...
         1.41, 2.40, 4.3888, 8.365 ...
         10.3522580, 40.245, 100.16, 1000.0276984];
   
sigma_normal_32 = [0.350, 0.360, 0.380, 0.420 ... 
         0.440, 0.540, 0.740, 1.139 ...
         1.335, 2.329677, 4.3198412, 8.30 ...
         10.2913492, 40.21079365, 100.13, 1000.02419047];

%------------------ Plot -------------------%
figure(1)
loglog(tau_shiba,sigma_shiba,"ko--")
hold on
loglog(tau,sigma_unpinned,"ks-")
hold on
loglog(tau,sigma_normal_8,"rs-")
hold on
loglog(tau,sigma_normal_16,"gs-")
hold on
loglog(tau,sigma_normal_32,"bs-")
hold on
grid on
xlabel("\tau (k_BT/a_0^2)")
ylabel("\sigma (k_BT/a_0^2)")
legend('Shiba et al.','Unpinned','8% pinning','16% pinning','32% pinning')

%------------------------------Refined---------------------------%

tau = [0.01, 0.08, 0.4, 1.0];

unpinned = [0.501      ,0.570      ,0.890     ,1.479];
normal_4 = [0.48492063 ,0.55451613 ,0.874     ,1.4659];
normal_8 = [0.4532     ,0.525      ,0.842     ,1.451];
normal_12 = [0.4503    ,0.519      ,0.8364    ,1.42913];
normal_16 = [0.431493  ,0.4994     ,0.814747  ,1.41275];
normal_20 = [0.410498  ,0.480      ,0.796215  ,1.39092];
normal_24 = [0.391599  ,0.461333   ,0.774648  ,1.37455];
normal_28 = [0.371088  ,0.44025629 ,0.7573912 ,1.3547773];
normal_32 = [0.349648  ,0.419573   ,0.741193  ,1.34209];

sig_1 = [0.501, 0.48188, 0.4532, 0.4503, 0.431493, 0.410498, 0.391599, 0.36664, 0.349648] / 0.501;
sig_2 = [0.570, 0.555, 0.525, 0.519, 0.4994, 0.480, 0.461333, 0.433778, 0.419573] / 0.570;
sig_3 = [0.890, 0.874, 0.842, 0.8364, 0.814747, 0.796215, 0.774648, 0.750836, 0.741193] / 0.890;
sig_4 = [1.479,1.4659,1.451,1.42913,1.41275,1.39092,1.37455,1.34782,1.34209] / 1.479;

sig_by_tau_1 = [0.501, 0.48188, 0.4532, 0.4503, 0.431493, 0.410498, 0.391599, 0.36664, 0.349648] / 0.01;
sig_by_tau_2 = [0.570, 0.555, 0.525, 0.519, 0.4994, 0.480, 0.461333, 0.433778, 0.419573] / 0.08;
sig_by_tau_3 = [0.890, 0.874, 0.842, 0.8364, 0.814747, 0.796215, 0.774648, 0.750836, 0.741193] / 0.4;
sig_by_tau_4 = [1.479,1.4659,1.451,1.42913,1.41275,1.39092,1.37455,1.34782,1.34209];

percentages = [0.0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32];

%------------------------------------------------------%
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
legend('Unpinned','4% pinning','8% pinning','12% pinning','16% pinning', ...
    '20% pinning','24% pinning','28% pinning','32% pinning', 'location','northwest');

%-------------------------------------------------------%
figure(3)
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
%-------------------------------------------------------%
figure(4)
plot(percentages, sig_by_tau_1,'Marker','o')
hold on
plot(percentages, sig_by_tau_2,'Marker','o')
hold on
plot(percentages, sig_by_tau_3,'Marker','o')
hold on
plot(percentages, sig_by_tau_4,'Marker','o')
hold on
grid on
xlabel(" Fraction  of pinning")
ylabel("\sigma/\tau")
grid on
legend('\tau=0.01 ','\tau=0.08','\tau=0.4','\tau=1.0');

%----------------------sigma Vs excess area------------------------%/
excess_area_unpinned = [0.0203,0.0202, 0.0194,0.0196,...
    0.0197,0.0194,0.0189,0.0181,...
    0.0170,0.0157,0.0133,0.0107,...
    0.0100,0.0056,...
    0.0032,...
    4.7*1e-4];
delta_a = 100*(excess_area_unpinned(1) - excess_area_unpinned);

excess_area_normal_8 = [0.0203,0.0202, 0.0194,0.0196,...
    0.0197,0.0194,0.0189,0.0181,...
    0.0170,0.0157,0.0133,0.0107,...
    0.0100,0.0056,...
    0.0032,...
    4.7*1e-4];
delta_a_normal_8 = 100*(excess_area_normal_8(1) - excess_area_normal_8);

excess_area_normal_16 = [0.0089,0.0089, 0.0089,0.0087,...
    0.0089,0.0088,0.0087,0.0086,...
    0.0085,0.0082,0.0078,0.0071,...
    0.0067,0.0043,...
    0.0025,...
    4.2367e-04];
delta_a_normal_16 = 100*(excess_area_normal_16(1) - excess_area_normal_16);

excess_area_normal_32 = [0.0059,0.0059, 0.0060,0.0059,...
    0.0060,0.0059,0.0059,0.0058,...
    0.0058,0.0057,0.0055,0.0052,...
    0.0049,0.0035,...
    7.4078e-04,...
    3.4295e-04];
delta_a_normal_32 = 100*(excess_area_normal_32(1) - excess_area_normal_32);


figure (5)
semilogy(delta_a, sigma_unpinned/sigma_unpinned(1),'bs-');
hold on
semilogy(delta_a_normal_8, sigma_normal_8/sigma_normal_8(1),'rs-');
hold on
semilogy(delta_a_normal_16, sigma_normal_16/sigma_normal_16(1),'gs-');
hold on
semilogy(delta_a_normal_32, sigma_normal_32/sigma_normal_32(1),'ks-');
xlabel(" \Deltaa = a_0-a (%)")
ylabel("\sigma/\sigma_0")
grid on
legend('Unpinned','8% pinning','16% pinning','32% pinning')
