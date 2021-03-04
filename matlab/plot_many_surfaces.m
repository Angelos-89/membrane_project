dir_pin_0  = "/home/angelos-89/runs/pinned/long_runs/pin_0/batch0";
dir_pin_10 = "/home/angelos-89/runs/pinned/long_runs/pin_10/batch0";
dir_pin_50 = "/home/angelos-89/runs/pinned/long_runs/pin_50/batch0";

hdata_pin_0  = h5read(dir_pin_0  + "/hfield_0.h5","/RectMesh")';
hdata_pin_10 = h5read(dir_pin_10 + "/hfield_0.h5","/RectMesh")';
hdata_pin_50 = h5read(dir_pin_50 + "/hfield_0.h5","/RectMesh")';

subplot(2,2,1)
surf(hdata_pin_0)
title('Unpinned membrane')
xlabel("x (a_0)")
ylabel("y (a_0)")
zlabel("z (a_0)")
zlim([-6,6])
grid on

subplot(2,2,2)
surf(hdata_pin_10)
title('Pinned membrane (10%)')
xlabel("x (a_0)")
ylabel("y (a_0)")
zlabel("z (a_0)")
zlim([-6,6])
grid on

subplot(2,2,[3,4])
surf(hdata_pin_50)
title('Pinned membrane (50%)')
xlabel("x (a_0)")
ylabel("y (a_0)")
zlabel("z (a_0)")
zlim([-6,6])
grid on