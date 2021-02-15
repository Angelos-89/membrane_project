%% import data
samples = LoadSamples("sampling0.txt");
height  = h5read("hfield0.h5","/RectMesh");
height  = inner(height,2);

%% lists of data
area   = samples(:,1);
alpha  = samples(:,3);
energy = samples(:,9);
 
%% plot data
% figure(1)
% surf(height)
% title("Membrane height field")
% xlabel("x");
% ylabel("y")
% zlabel("z (a_0)")
% zlim([-6,6])

figure(2)
plot(energy);
title("Total energy vs MC steps");
xlabel("MC steps");
ylabel("Energy (k_bT/N^2)");

% figure(3)
% plot(area);
% title("Total area vs MC steps");
% xlabel("MC steps");
% ylabel("Area (a_0^2/N^2)");
% 
% figure(4)
% plot(alpha);
% title("Lattice spacing vs MC steps");
% xlabel("MC steps");
% ylabel("alpha (a_0)");
% 
% %% calculate means and stds and plot histograms
% begin = 2*1e4;
% 
% area = area(begin:end);
% area_mean = mean(area);
% area_std = std(area);
% area_pd = normpdf(area,area_mean,area_std)/100;
% figure(5)
% histfit(area)
% title("Area per degrees of freedom");
% xlabel("Area per degrees of freedom (a_0^2/N^2)");
% ylabel("Occurances");
% grid on
% 
% alpha = alpha(begin:end);
% alpha_mean = mean(alpha);
% alpha_std = std(alpha);
% alpha_pd = normpdf(alpha,alpha_mean,alpha_std)/100;
% figure(6)
% histfit(alpha);
% title("Lattice spacing");
% xlabel("Lattice spacing (a_0)");
% ylabel("Occurances");
% grid on
% 
% energy = energy(begin:end);
% energy_mean = mean(energy);
% energy_std = std(energy);
% energy_pd = normpdf(energy,energy_mean,energy_std)/100;
% figure(7)
% histfit(energy)
% title("Energy per degrees of freedom");
% xlabel("Energy per degree of freedom (k_bT/N^2)");
% ylabel("Occurances");
% grid on

%% define functions
function samples = LoadSamples(filename)
samples = importdata(filename); 
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