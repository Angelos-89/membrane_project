tau = [0.01, 0.1, 1.0, 10., 100., 1000.];
exc_area_pin_0  = [0.0677, 0.0656, 0.0488, 0.0216, 0.0044, 4.67*1e-4];
exc_area_pin_10 = [0.0403, 0.0405, 0.0346, 0.0183, 0.0038, 4.48*1e-4];
exc_area_pin_50 = [0.0145, 0.0144, 0.0138, 0.0090, 0.0022, 2.95*1e-4];

figure(1)
semilogy(exc_area_pin_0,tau,'r*-');
hold on
semilogy(exc_area_pin_10,tau,'bo-');
hold on
semilogy(exc_area_pin_50,tau,'go-');
title("Frame tension \tau Vs excess area");
xlabel("\langle{A-A_p}/{A}\rangle");
ylabel("\tau (k_B/a_0^2)");
legend('unpinned','10% pinning','50% pinning','location','northwest');
grid on;
