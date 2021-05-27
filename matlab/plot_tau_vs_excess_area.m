DATA = importdata("excess-area-data.txt");
tau = DATA.data(:,1);
ExcessAreas_unpinned = DATA.data(:,2);
ExcessAreas_8 = DATA.data(:,4);
ExcessAreas_16 = DATA.data(:,6);
ExcessAreas_32 = DATA.data(:,8);

%----------------------
figure(1)
semilogy(ExcessAreas_unpinned*100,tau,"ro-");
hold on
semilogy(ExcessAreas_8*100,tau,"bo-");
hold on
semilogy(ExcessAreas_16*100,tau,"go-");
hold on
semilogy(ExcessAreas_32*100,tau,"ko-");
ylabel("$\tau (\frac{k_BT}{a_0^2})$",'interpreter','latex','FontSize',12);
xlabel("Excess area $\%$",'interpreter','latex','FontSize',12);
grid on
legend('unpinned','8% pinning','16% pinning','32% pinning')
%----------------------

sigmas_unpinned = [0.506143, 0.515, 0.5348571, 0.5728571, 0.5938709, 0.69, 0.89, 1.2805 ...
                   1.479, 2.465, 4.448, 8.41, 10.4, 40.28, 100.18, 1000.03];

sigmas_8 = [0.46928,0.47666, 0.4972, 0.5359677, 0.5594, 0.653, 0.85666, 1.2496825...
            1.451, 2.4396825, 4.42, 8.39066, 10.379677, 40.264, 100.17, 1000.03];
        
sigmas_16 = [0.427, 0.438, 0.4602, 0.496, 0.519333, 0.615, 0.816, 1.21333, 1.41, 2.400 ...
            4.388, 8.365, 10.3522580, 40.245, 100.16, 1000.0276984];
               
sigmas_32 = [0.350,0.36,0.38,0.42,0.440,0.540,0.74,1.139,1.335, 2.329677, 4.3198412, 8.3 ...
            10.2913492, 40.21079365, 100.13, 1000.02419047];
               
sByS0_unpinned = sigmas_unpinned/sigmas_unpinned(1);
sByS0_8 = sigmas_8/sigmas_unpinned(1);
sByS0_16 = sigmas_16/sigmas_unpinned(1);
sByS0_32 = sigmas_32/sigmas_unpinned(1);

a0_unpinned = ExcessAreas_unpinned(1);
a0_8 = ExcessAreas_8(1);
a0_16 = ExcessAreas_16(1);
a0_32 = ExcessAreas_32(1);

aa_unpinned = a0_unpinned - ExcessAreas_unpinned;
aa_8 = a0_8 - ExcessAreas_8;
aa_16 = a0_16 - ExcessAreas_16;
aa_32 = a0_32 - ExcessAreas_32;

figure(2)
semilogy(aa_unpinned*100,sByS0_unpinned,'ro-');
hold on
semilogy(aa_8*100,sByS0_unpinned,'bo-');
hold on
semilogy(aa_16*100,sByS0_unpinned,'go-');
hold on
semilogy(aa_32*100,sByS0_unpinned,'ko-');
ylabel("$\sigma (\frac{k_BT}{a_0^2})$",'interpreter','latex','FontSize',12);
xlabel("$\alpha_0 - \alpha \%$",'interpreter','latex','FontSize',12);
grid on
legend('unpinned','8% pinning','16% pinning','32% pinning')
