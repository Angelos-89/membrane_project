tau = [0.01, 0.1, 1.0, 10., 100., 1000.];
sigma_pin_0 = [0.500,0.570,1.440,10.280,100.070,1000.008];
sigma_pin_10 = [0.453,0.5353,1.408,10.263,100.060,1000.01];
sigma_pin_50 = [0.264,0.350,1.30,10.16,100.040,1000.0057];

exc_area_pin_0  = [0.0757, 0.0711, 0.0488, 0.0216, 0.0044, 4.67*1e-4];
exc_area_pin_10 = [0.0403, 0.0405, 0.0346, 0.0183, 0.0038, 4.48*1e-4];
exc_area_pin_50 = [0.0145, 0.0144, 0.0138, 0.0090, 0.0022, 2.95*1e-4];

figure(1)
semilogy(exc_area_pin_0,tau,'ro-');
hold on
semilogy(exc_area_pin_0,sigma_pin_0,'ro--','MarkerFaceColor', 'r','HandleVisibility','off');
hold on
semilogy(exc_area_pin_10,tau,'bo-');
hold on
semilogy(exc_area_pin_10,sigma_pin_10,'bo--','MarkerFaceColor', 'b','HandleVisibility','off');
hold on
semilogy(exc_area_pin_50,tau,'go-');
hold on
semilogy(exc_area_pin_50,sigma_pin_50,'go--','MarkerFaceColor', 'g','HandleVisibility','off');
xlabel("\langle{A-A_p}/{A}\rangle");
ylabel("\tau,\sigma(k_B/a_0^2)");
legend('unpinned','10% pinning','50% pinning','location','northwest');
grid on;
