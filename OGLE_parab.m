data = load("~/Desktop/phot.dat");
x = data(:,1);
magnitude = data(:,2);
flux = 1*10.^(0.4*(-magnitude));
disp(flux)
err = data(:,3);
sigma_flux = -0.4*log(10).*err.*10.^(0.4*(-magnitude));

flux_val = [];
x_val = [];
sigma_val = [];

for i = 1:60
    flux_val(i)= flux(i+470);
    x_val(i)= x(i+470);
    sigma_val(i)= sigma_flux(i+470);
end


[a,b,c] = parab(x_val,flux_val,sigma_val);
%[a,b,c] = parab(x,flux,sigma_flux);
