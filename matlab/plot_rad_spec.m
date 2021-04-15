load("hfield_spec_0.txt");
k  = hfield_spec_0(2:40,1);
sp = hfield_spec_0(2:40,2);
sp = sp*6400;
% figure(1)
% plot(k,sp,"bo-")
% grid on
figure(2)
plot(log(k(1:end)),log(sp(1:end)),"bo-");
hold on
f=fit(log(k(1:end)),log(sp(1:end)),'poly1');
plot(f,log(k),log(sp),"b");
hold on
plot(log(k),-log(0.5/2/pi*k.^3),"g");

hold on
plot(log(k),-log(0.5/2/pi*k),"k--");


% hold on
% g=fit(log(k(4:end)),log(1./sp(4:end)),'poly1');
% plot(g,log(k),log(1./sp),"k--");
% hold on
% plot(log(k),log(10/2/pi*k.^3),"g--");
% xlim([-1,1])
grid on
