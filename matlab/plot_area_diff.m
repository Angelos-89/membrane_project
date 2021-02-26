u_dir = "/home/angelos-89/runs/pinned/long_runs/pin_0/tau_0.01";
p_dir = "/home/angelos-89/runs/pinned/long_runs/pin_10/tau_0.01";
filename = "timeseries_0.txt";
plot_area_difference(u_dir,p_dir,filename);

function plot_area_difference(u_dir, p_dir, filename)
    u_data = load(u_dir + "/" + filename);
    p_data = load(p_dir + "/" + filename);
    u_area = u_data(:,3);
    u_prj = u_data(:,4);
    plot( (u_area-u_prj)./u_area,"b");
    hold on
    p_area = p_data(:,3);
    p_prj = p_data(:,4);
    plot( (p_area-p_prj)./p_area,"r");
    legend("unpinned", "pinned");
    xlabel("MC steps");
    ylabel("(A-A_p) / A")
    grid on
end