DATA = importdata("excess-area.data");
tau = DATA.data(:,1);

sigmas_unpinned = [0.506143, 0.515, 0.5348571, 0.5728571, 0.5938709, 0.69, 0.89, 1.2805 ...
                   1.479, 2.465, 4.448, 8.41, 10.4, 40.28, 100.18, 1000.03];

sigmas_8 = [0.46928,0.47666, 0.4972, 0.5359677, 0.5594, 0.653, 0.85666, 1.2496825...
            1.451, 2.4396825, 4.42, 8.39066, 10.379677, 40.264, 100.17, 1000.03];
        
sigmas_16 = [0.427, 0.438, 0.4602, 0.496, 0.519333, 0.615, 0.816, 1.21333, 1.41, 2.400 ...
            4.388, 8.365, 10.3522580, 40.245, 100.16, 1000.0276984];
               
sigmas_32 = [0.350,0.36,0.38,0.42,0.440,0.540,0.74,1.139,1.335, 2.329677, 4.3198412, 8.3 ...
            10.2913492, 40.21079365, 100.13, 1000.02419047];
               
ExcessAreas_unpinned = DATA.data(:,2);
ExcessAreas_8 = DATA.data(:,5);
ExcessAreas_16 = DATA.data(:,8);
ExcessAreas_32 = DATA.data(:,11);
 
error_bar_min_unpinned = (ExcessAreas_unpinned - DATA.data(:,3));
error_bar_min_8 = (ExcessAreas_8 - DATA.data(:,6));
error_bar_min_16 = (ExcessAreas_16 - DATA.data(:,9));
error_bar_min_32 = (ExcessAreas_32 - DATA.data(:,12));

error_bar_max_unpinned = (DATA.data(:,4) - ExcessAreas_unpinned);
error_bar_max_8 = (DATA.data(:,7) - ExcessAreas_8);
error_bar_max_16 = (DATA.data(:,10) - ExcessAreas_16);
error_bar_max_32 = (DATA.data(:,13) - ExcessAreas_32);

%--------------------------------------------------------------------------
figure(1)
errorbar(ExcessAreas_unpinned*100 , tau , error_bar_min_unpinned*100 , error_bar_max_unpinned*100 , 'horizontal' , "ro-");
hold on
errorbar(ExcessAreas_8*100 , tau , error_bar_min_8*100 , error_bar_max_8*100 , 'horizontal' , 'bo-');
hold on
errorbar(ExcessAreas_16*100 , tau , error_bar_min_16*100 , error_bar_max_16*100 , 'horizontal' , 'go-');
hold on
errorbar(ExcessAreas_32*100 , tau,error_bar_min_32*100 , error_bar_max_32*100 , 'horizontal' , 'ko-');
ylabel("$\tau \ (\frac{k_BT}{a_0^2})$",'interpreter','latex','FontSize',14);
xlabel("Excess area $\frac{\langle A \rangle - \langle A_p \rangle}{\langle A \rangle} \ \%$",'interpreter','latex','FontSize',13);
grid on
legend('unpinned','8% pinning','16% pinning','32% pinning')
set(gca,'YScale','log');
%--------------------------------------------------------------------------

sByS0_unpinned = sigmas_unpinned/sigmas_unpinned(1);
sByS0_8 = sigmas_8/sigmas_8(1);
sByS0_16 = sigmas_16/sigmas_16(1);
sByS0_32 = sigmas_32/sigmas_32(1);

a0_unpinned = ExcessAreas_unpinned(1);
a0_8 = ExcessAreas_8(1);
a0_16 = ExcessAreas_16(1);
a0_32 = ExcessAreas_32(1);

aa_unpinned = (a0_unpinned - ExcessAreas_unpinned);
aa_8 = (a0_8 - ExcessAreas_8);
aa_16 = (a0_16 - ExcessAreas_16);
aa_32 = (a0_32 - ExcessAreas_32);

%error bars
ExcessAreas_unpinned_min = DATA.data(:,3);
ExcessAreas_8_min = DATA.data(:,6);
ExcessAreas_16_min = DATA.data(:,9);
ExcessAreas_32_min = DATA.data(:,12);

ExcessAreas_unpinned_max = DATA.data(:,4);
ExcessAreas_8_max = DATA.data(:,7);
ExcessAreas_16_max = DATA.data(:,10);
ExcessAreas_32_max = DATA.data(:,13);

a0_unpinned_min = DATA.data(1,3);
a0_8_min = DATA.data(1,6);
a0_16_min = DATA.data(1,9);
a0_32_min = DATA.data(1,12);

a0_unpinned_max = DATA.data(1,4);
a0_8_max = DATA.data(1,7);
a0_16_max = DATA.data(1,10);
a0_32_max = DATA.data(1,13);

aa_unpinned_min = (a0_unpinned_min - ExcessAreas_unpinned_min);
aa_8_min = (a0_8_min - ExcessAreas_8_min);
aa_16_min = (a0_16_min - ExcessAreas_16_min);
aa_32_min = (a0_32_min - ExcessAreas_32_min);

aa_unpinned_max = (a0_unpinned_max - ExcessAreas_unpinned_max);
aa_8_max = (a0_8_max - ExcessAreas_8_max);
aa_16_max = (a0_16_max - ExcessAreas_16_max);
aa_32_max = (a0_32_max - ExcessAreas_32_max);

error_bars_unpinned_min = abs(aa_unpinned - aa_unpinned_min);  
error_bars_8_min = abs(aa_8 - aa_8_min);
error_bars_16_min = abs(aa_16 - aa_16_min);
error_bars_32_min = abs(aa_32 - aa_32_min);

error_bars_unpinned_max = abs(aa_unpinned_max - aa_unpinned);  
error_bars_8_max = abs(aa_8_max - aa_8);
error_bars_16_max = abs(aa_16_max - aa_16);
error_bars_32_max = abs(aa_32_max - aa_32);

%---------------------------------------------%
figure(2)
errorbar(aa_unpinned*100 , sByS0_unpinned,error_bars_unpinned_min*100 , error_bars_unpinned_max*100 , 'horizontal','ro');
hold on
errorbar(aa_8*100 , sByS0_8 , error_bars_8_min*100 , error_bars_8_max*100,'horizontal','bo');
hold on
errorbar(aa_16*100 , sByS0_16 , error_bars_16_min*100 , error_bars_16_max*100 ,'horizontal','go');
hold on
errorbar(aa_32*100 , sByS0_32 , error_bars_32_min*100 , error_bars_32_max*100 ,'horizontal','ko');
grid on
ylabel("$\frac{\sigma}{\sigma_0} (\frac{k_BT}{a_0^2})$",'interpreter','latex','FontSize',16);
xlabel("$ (\alpha_0 - \alpha)\ \% $",'interpreter','latex','FontSize',12);
set(gca,'YScale','log');
legend('unpinned','8% pinning','16% pinning','32% pinning')
xlim([0,3.0])

a2 = axes(figure(2));
a2.Position = [0.58 0.15 0.3 0.23]; % xlocation, ylocation, xsize, ysize
errorbar(aa_unpinned(1:6)*100,sByS0_unpinned(1:6),error_bars_unpinned_min(1:6)*100,error_bars_unpinned_max(1:6)*100,'horizontal','ro'); axis tight;
hold on
errorbar(aa_8(1:6)*100,sByS0_8(1:6),error_bars_8_min(1:6)*100,error_bars_8_max(1:6)*100,'horizontal','bo'); axis tight;
hold on
errorbar(aa_16(1:6)*100,sByS0_16(1:6),error_bars_16_min(1:6)*100,error_bars_16_max(1:6)*100,'horizontal','go'); axis tight;
hold on
errorbar(aa_32(1:6)*100,sByS0_32(1:6),error_bars_32_min(1:6)*100,error_bars_32_max(1:6)*100,'horizontal','ko'); axis tight;
grid on;
% annotation('arrow',[.52 .2],[.16 .12])
grid on