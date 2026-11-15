data = load("~/Desktop/phot.dat");
x = data(:,1);
magnitude = data(:,2);
flux = 1*10.^(0.4*(-magnitude));
%disp(flux);
err = data(:,3);
sigma_flux = -0.4*log(10).*err.*10.^(0.4*(-magnitude));
num = linspace(1,574,574);
disp(x(514))
plot(num,flux);
xlabel('time');
ylabel('flux');

%MAX in i = 514
