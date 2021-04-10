%% import data
datastruct = importdata("timeseries_0.txt");
samples = LoadSamples(datastruct);
height  = h5read("hfield_0.h5","/RectMesh");
height  = inner(height,2);
begin = 2e3;

%% lists of data
Area     = samples(:,3);
Prj_area = samples(:,4);
Alpha    = samples(:,5);
Energy   = samples(:,9);
 
%% plot data

%-----------------------height--------------------------%
hdata      = height(:) - mean(height(:));
hdata_mean = mean(hdata);
hdata_std  = std(hdata);
hdata_pd   = normpdf(hdata,hdata_mean,hdata_std)/100;

figure(5)
histogram(hdata,'DisplayStyle','bar');
% title("Height field values");
xlabel("Height (a_0)");
ylabel("Occurances");
grid on
hold on
xline(hdata_mean,'r')

figure(1)
surf(height)
% title("Membrane height field")
xlabel("x (a_0)");
ylabel("y (a_0)")
zlabel("z (a_0)")
zlim([-7,7])
grid on

%------------------------energy-------------------------%
energy      = Energy(begin:end);
energy_mean = mean(energy);
energy_std  = std(energy);
energy_pd   = normpdf(energy,energy_mean,energy_std)/100;

figure(9)
histogram(energy,'DisplayStyle','stairs')
% title("Energy per degrees of freedom");
xlabel("Energy per degrees of freedom (k_BT/N^2)");
ylabel("Occurances");
grid on
hold on
xline(energy_mean,'r')

figure(2)
plot(Energy);
% title("Total energy vs MC steps");
xlabel("MC time");
ylabel("Energy (k_BT/N^2)");
grid on
hold on
yline(energy_mean,'r');

%------------------------area---------------------------%
area      = Area(begin:end);
area_mean = mean(area);
area_std  = std(area);
area_pd   = normpdf(area,area_mean,area_std)/100;

figure(6)
histogram(area,'DisplayStyle','stairs');
% title("Area per degrees of freedom");
xlabel("Area per degrees of freedom (a_0^2/N^2)");
ylabel("Occurances");
grid on
hold on
xline(area_mean,'r')

%------------area---------------%
figure(3)
plot(Area);
% title("Total area vs MC steps");
xlabel("MC time");
ylabel("Area (a_0^2/N^2)");
grid on
hold on
yline(area_mean,'r')

%------------------lattice spacing----------------------%
alpha      = Alpha(begin:end);
alpha_mean = mean(alpha);
alpha_std  = std(alpha);
alpha_pd   = normpdf(alpha,alpha_mean,alpha_std)/100;

figure(8)
histogram(alpha,'DisplayStyle','stairs');
% title("Lattice spacing");
xlabel("Lattice spacing (a_0)");
ylabel("Occurances");
grid on
hold on
xline(alpha_mean,'r')

figure(4)
plot(Alpha);
% title("Lattice spacing vs MC steps");
xlabel("MC time");
ylabel("alpha (a_0)");
ylim([alpha_mean-10*alpha_std,alpha_mean+10*alpha_std]);
grid on
hold on
yline(alpha_mean,'r')

%---------------------prj_area--------------------------%
prj_area      = Prj_area(begin:end);
prj_area_mean = mean(prj_area);
prj_area_std  = std(prj_area);
prj_area_pd   = normpdf(prj_area,prj_area_mean,prj_area_std)/100;

%--------------------excess area------------------------%
exc_area = (Area-Prj_area)./Area;
exc_area_mean = mean(exc_area);
exc_area_std  = std(exc_area);
exc_area_pd   = normpdf(exc_area,exc_area_mean,exc_area_std)/100;

figure(7)
plot(exc_area);
% title("Excess area vs MC steps");
xlabel("MC time");
ylabel("Excess area \langle{A-A_p}/{A}\rangle");
grid on
hold on
yline(exc_area_mean,'r');

%--------------------display excess area mean-----------%
disp((area_mean-prj_area_mean)/area_mean);

%% define functions
function samples = LoadSamples(struct)
    samples = struct.data; 
%     if length(samples(:,1)) < 1e4
%         error("Samples must be more than 1e4. Exiting.")
%     end
end

function field = inner(hfield,nghost)
    field = hfield'; 
    field(1:nghost,:)         = [];
    field(end-nghost+1:end,:) = [];
    field(:,1:nghost)         = [];
    field(:,end-nghost+1:end) = [];
end