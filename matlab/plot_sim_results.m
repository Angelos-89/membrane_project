%% import data
datastruct = importdata("timeseries_0.txt");
samples = datastruct.data;
height = h5read("hfield_0.h5","/RectMesh");
height = inner(height,2);
begin = 2e3;

%% lists of data
Area     = samples(:,3);
Prj_area = samples(:,4);
Alpha    = samples(:,5);
Energy   = samples(:,9);
 
%% plot data
%--------------------------height-----------------------------%
hdata = height(:) - mean(height(:));
hdata_mean = mean(hdata);
hdata_std  = std(hdata);

Xh_step = 1e-2;
Xh = hdata_mean-5*hdata_std : Xh_step : hdata_mean+5*hdata_std; 
hdata_pd = normpdf(Xh, hdata_mean, hdata_std) * Xh_step;

figure(1)
semilogy(Xh,hdata_pd,'ko')
xlabel("Height ($a_0$)",'interpreter','latex','FontSize',15);
ylabel("Normalized PDF",'interpreter','latex','FontSize',15);
xline(hdata_mean,'g')
xline(hdata_mean - hdata_std,'r--')
xline(hdata_mean + hdata_std,'r--')
xline(hdata_mean - 2*hdata_std,'b--')
xline(hdata_mean + 2*hdata_std,'b--')
legend('PDF','mean','std')
grid on

figure(2)
surf(height)
xlabel("x ($a_0$)",'interpreter','latex','FontSize',15);
ylabel("y ($a_0$)",'interpreter','latex','FontSize',15);
zlabel("z ($a_0$)",'interpreter','latex','FontSize',15);
zlim([-7,7])
grid on
 
%------------------------energy-------------------------%
energy = Energy(begin:end);
energy_mean = mean(energy);
energy_std = std(energy);

Xe_step = 1e-3;
Xe = energy_mean-5*energy_std : Xe_step : energy_mean+5*energy_std; 
energy_pd = normpdf(Xe, energy_mean, energy_std) * Xe_step;

figure(3)
loglog(Xe,energy_pd,'ko')
xlabel("Energy per degrees of freedom ($\frac{k_BT}{N^2}$)",'interpreter','latex','FontSize',15);
ylabel("Normalized PDF",'interpreter','latex','FontSize',15);
xline(energy_mean,'g')
xline(energy_mean - energy_std,'r--')
xline(energy_mean + energy_std,'r--')
xline(energy_mean - 2*energy_std,'b--')
xline(energy_mean + 2*energy_std,'b--')
legend('PDF','mean','std')
grid on
 
figure(4)
plot(Energy);
xlabel("Samples",'interpreter','latex','FontSize',15);
ylabel("Energy ($\frac{k_BT}{N^2}$)",'interpreter','latex','FontSize',15);
grid on
hold on
yline(energy_mean,'r');
  
%------------------------area----------------------------%
area = Area(begin:end);
area_mean = mean(area);
area_std = std(area);
 
Xa_step = 1e-3;
Xa = area_mean-5*area_std : Xa_step : area_mean+5*area_std; 
area_pd = normpdf(Xa, area_mean, area_std) * Xa_step;
 
figure(5)
loglog(Xa,area_pd,'ko')
xlabel("Area per degrees of freedom ($\frac{a_0^2}{N^2}$)",'interpreter','latex','FontSize',15);
ylabel("Normalized PDF",'interpreter','latex','FontSize',15);
xline(area_mean,'g')
xline(area_mean - area_std,'r--')
xline(area_mean + area_std,'r--')
xline(area_mean - 2*area_std,'b--')
xline(area_mean + 2*area_std,'b--')
legend('PDF','mean','std')
grid on

figure(6)
plot(Area);
xlabel("Samples",'interpreter','latex','FontSize',15);
ylabel("Area ($\frac{a_0^2}{N^2}$)",'interpreter','latex','FontSize',15);
grid on
hold on
yline(area_mean,'r')

%------------------lattice spacing----------------------%
alpha      = Alpha(begin:end);
alpha_mean = mean(alpha);
alpha_std  = std(alpha);

Xalph_step = 1e-3;
Xalph = alpha_mean-5*alpha_std : Xalph_step : alpha_mean+5*alpha_std; 
alpha_pd = normpdf(Xalph, alpha_mean, alpha_std) * Xalph_step;
 
figure(7)
loglog(Xalph,alpha_pd,'ko')
xline(alpha_mean,'g')
xline(alpha_mean - alpha_std,'r--')
xline(alpha_mean + alpha_std,'r--')
xline(alpha_mean - 2*alpha_std,'b--')
xline(alpha_mean + 2*alpha_std,'b--')
legend('PDF','mean','std')
grid on
xlabel("Lattice spacing ($a_0$)",'interpreter','latex','FontSize',15);
ylabel("Normalized PDF",'interpreter','latex','FontSize',15);
grid on
hold on
 
figure(8)
plot(Alpha);
xlabel("Samples",'interpreter','latex','FontSize',15);
ylabel("Lattice spacing ($a_0$)",'interpreter','latex','FontSize',15);
ylim([alpha_mean-10*alpha_std,alpha_mean+10*alpha_std]);
grid on
hold on
yline(alpha_mean,'r')
 
%---------------------prj_area--------------------------%
prj_area = Prj_area(begin:end);
prj_area_mean = mean(prj_area);
prj_area_std = std(prj_area);

Xprj_step = 1e-3;
Xprj = prj_area_mean-5*prj_area_std : Xprj_step : prj_area_mean+5*prj_area_std; 
prj_area_pd = normpdf(Xprj, prj_area_mean, prj_area_std) * Xprj_step;

figure(9)
loglog(Xprj,prj_area_pd,'ko')
xline(prj_area_mean,'g')
xline(prj_area_mean - prj_area_std,'r--')
xline(prj_area_mean + prj_area_std,'r--')
xline(prj_area_mean - 2*prj_area_std,'b--')
xline(prj_area_mean + 2*prj_area_std,'b--')
legend('PDF','mean','std')
grid on
xlabel("Projected area ($a_0^2$)",'interpreter','latex','FontSize',15);
ylabel("Normalized PDF",'interpreter','latex','FontSize',15);
grid on
hold on
 
figure(10)
plot(prj_area);
xlabel("Samples",'interpreter','latex','FontSize',15);
ylabel("Projected area ($a_0^2$)",'interpreter','latex','FontSize',15);
grid on
hold on
yline(alpha_mean,'r')
 
%--------------------excess area------------------------%
exc_area = (Area-Prj_area)./Area;
exc_area_mean = mean(exc_area);
exc_area_std = std(exc_area);
 
figure(11)
plot(exc_area);
xlabel("Samples",'interpreter','latex','FontSize',15);
ylabel("Excess area $\frac{A-A_p}{A}$",'interpreter','latex','FontSize',15);
% set(hh, 'FontSize', 15) 
% set(hh,'FontWeight','bold') %bold font
grid on
hold on
yline(exc_area_mean,'r');
 
%--------------------display excess area mean-----------%
disp("Excess area is: " + (area_mean-prj_area_mean)/area_mean);
 
%% define functions
function field = inner(hfield,nghost)
    field = hfield'; 
    field(1:nghost,:)         = [];
    field(end-nghost+1:end,:) = [];
    field(:,1:nghost)         = [];
    field(:,end-nghost+1:end) = [];
end
