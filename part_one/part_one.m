data = load("~/Desktop/phot.dat");
time = data(:,1);
magnitude = data(:,2);
% Converting magnitude to flux based on eq.(1.14) when F0=1 and m0=0
flux = 1*10.^(0.4*(-magnitude));
err = data(:,3);
% Converting magnitude to flux based on the derivative of eq.(1.14) with respect to m, when F0=1 and m0=0
sigma_flux = -0.4*log(10).*err.*10.^(0.4*(-magnitude));

% Normalizing the flux values with respect to flux value without
% microlensing.
flux = flux/(5*10^-9);
sigma_flux = sigma_flux/(5*10^-9);

% Slicing a small interval out of OGLE's data that behaves approximately
% like a parabola.

flux_val = [];
time_val = [];
sigma_val = [];

for i = 1:40
    flux_val(i) = flux(i+493);
    time_val(i) = time(i+493);
    sigma_val(i) = sigma_flux(i+493);
end

% Normalizing the values of time with respect to the first value.
time_norm_factor = time_val(1);
time_val = time_val-time_norm_factor;

% Calculating the values of a,b,c using parab function.
[a,b,c] = parab(time_val,flux_val,sigma_val);

parabola = a.*time_val.^2+b.*time_val+c;

% Time of measuring the maximum value of flux in Julian day units.
min_time = -b/(2*a)+time_norm_factor;
disp(min_time);

flux_max = a*(min_time-time_norm_factor)^2+b*(min_time-time_norm_factor)+c;
disp(flux_max);

% Plot
plot(time_val,parabola);
hold on
errorbar(time_val,flux_val,sigma_val, ".");
hold off
