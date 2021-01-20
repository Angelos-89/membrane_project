%% import data
height = h5read("hfield.h5","/RectMesh");
% height2 = h5read("height_fast_way.h5","/Mesh2D");

% hdiff = abs(height-height2); 

sampling = importdata("sampling.txt");


energy = sampling(:,8);
area = sampling(:,1); 
alpha = sampling(:,3);
 
% %% plot data

figure(1)
surf(height)
title("Membrane height field")
xlabel("x");
ylabel("y")
zlabel("z (a_0)")
zlim([-6,6])

figure(2)
plot(energy);
title("Total energy vs MC steps");
xlabel("MC steps");
ylabel("Energy (k_bT/N^2)");

figure(3)
plot(area);
title("Total area vs MC steps");
xlabel("MC steps");
ylabel("Area (a_0^2/N^2)");

figure(4)
plot(alpha);
title("Lattice spacing vs MC steps");
xlabel("MC steps");
ylabel("alpha (a_0)");
