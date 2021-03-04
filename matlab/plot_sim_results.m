%% import data
datastruct = importdata("timeseries_0.txt");
samples = LoadSamples(datastruct);
height  = h5read("hfield_0.h5","/RectMesh");
height  = inner(height,2);

%% lists of data
area     = samples(:,3);
prj_area = samples(:,4);
alpha    = samples(:,5);
energy   = samples(:,9);
 
%% plot data

%------------height-------------%
figure(1)
surf(height)
title("Membrane height field")
xlabel("x");
ylabel("y")
zlabel("z (a_0)")
zlim([-6,6])
grid on

% %------------energy-------------%
figure(2)
plot(energy);
title("Total energy vs MC steps");
xlabel("MC steps");
ylabel("Energy (k_bT/N^2)");
grid on

%------------area---------------%
figure(3)
plot(area);
title("Total area vs MC steps");
xlabel("MC steps");
ylabel("Area (a_0^2/N^2)");
grid on

%------------alpha--------------%
figure(4)
plot(alpha);
title("Lattice spacing vs MC steps");
xlabel("MC steps");
ylabel("alpha (a_0)");
grid on

%% calculate means, stds and plot histograms
begin = 1e5;

%-----------------------height--------------------------%
hdata      = height(:) - mean(height(:));
hdata_mean = mean(hdata);
hdata_std  = std(hdata);
hdata_pd   = normpdf(hdata,hdata_mean,hdata_std)/100;

figure(5)
histogram(hdata,'DisplayStyle','bar');
title("Height field values");
xlabel("Height (a_0)");
ylabel("Occurances");
grid on

%------------------------area---------------------------%
area      = area(begin:end);
area_mean = mean(area);
area_std  = std(area);
area_pd   = normpdf(area,area_mean,area_std)/100;

figure(6)
histogram(area,'DisplayStyle','stairs');
title("Area per degrees of freedom");
xlabel("Area per degrees of freedom (a_0^2/N^2)");
ylabel("Occurances");
grid on

%---------------------prj_area--------------------------%
prj_area      = prj_area(begin:end);
prj_area_mean = mean(prj_area);
prj_area_std  = std(prj_area);
prj_area_pd   = normpdf(prj_area,prj_area_mean,prj_area_std)/100;

%--------------------excess area------------------------%
exc_area = (area-prj_area)./area;
figure(7)
plot(exc_area);
title("Excess area vs MC steps");
xlabel("MC steps");
ylabel("Excess area (A-A_p) / A");
grid on

exc_area_mean = mean(exc_area);
exc_area_std  = std(exc_area);
exc_area_pd   = normpdf(exc_area,exc_area_mean,exc_area_std)/100;

%------------------lattice spacing----------------------%
alpha      = alpha(begin:end);
alpha_mean = mean(alpha);
alpha_std  = std(alpha);
alpha_pd   = normpdf(alpha,alpha_mean,alpha_std)/100;

figure(8)
histogram(alpha,'DisplayStyle','stairs');
title("Lattice spacing");
xlabel("Lattice spacing (a_0)");
ylabel("Occurances");
grid on

%------------------------energy-------------------------%
energy      = energy(begin:end);
energy_mean = mean(energy);
energy_std  = std(energy);
energy_pd   = normpdf(energy,energy_mean,energy_std)/100;
figure(9)
histogram(energy,'DisplayStyle','stairs')
title("Energy per degrees of freedom");
xlabel("Energy per degrees of freedom (k_bT/N^2)");
ylabel("Occurances");
grid on

%% define functions
function samples = LoadSamples(struct)
    samples = struct.data; 
    if length(samples(:,1)) < 1e4
        error("Samples must be more than 1e4. Exiting.")
    end
end

function field = inner(hfield,nghost)
    field = hfield'; 
    field(1:nghost,:)         = [];
    field(end-nghost+1:end,:) = [];
    field(:,1:nghost)         = [];
    field(:,end-nghost+1:end) = [];
end