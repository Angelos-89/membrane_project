load("hfield_spec_0.txt");
k  = hfield_spec_0(1:end-1,1);
sp = hfield_spec_0(1:end-1,2);
figure(1)
plot(k,sp,"bo-")
grid on
figure(2)
plot(log(k),log(sp),"bo-");
grid on
