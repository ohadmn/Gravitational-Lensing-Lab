% Hardcoded data per incident.
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
f_norm = flux(1);
flux = flux/f_norm;
sigma_flux = sigma_flux/(f_norm);
u_min = 0.121;
tau = 51.032;
T0 = 2.454476305378268e+06;

num_measurements = 1000;

% Defining the spectrum of u_min's values
u_min_lower_limit = u_min*0.5;
delta_u_min = (u_min*0.5)/(num_measurements/2);
% Defining the spectrum of tau's values
tau_lower_limit = tau*0.5;
delta_tau = (tau*0.5)/(num_measurements/2);

% The dimensions of chi_mat are the number of measurements.
chi_mat = zeros(num_measurements,num_measurements);

% Calculating chi squared with different variations from tau and u_min
% spectrum.
for i = 1:num_measurements
    for j = 1:num_measurements           
        chi_squared = 0;
        for k = 1:574
            t = (time(k)-T0)/(tau_lower_limit + delta_tau * i);
            u = ((u_min_lower_limit + delta_u_min * j) ^ 2 + t ^ 2) ^ (1/2);
            mu = (2 + u ^ 2) / (u * (4 + u ^ 2) ^ (1/2));
            
            chi_mat(i,j) = chi_mat(i,j) + (((flux(k) - mu) / sigma_flux(k)) ^ 2);
        end
    end
end
%disp(chi_mat);

% Finding the minimazation of the matrix.
[x, y] = matrix_minimization(chi_mat);
disp(x);
disp(y);

tau_final = tau_lower_limit+x*delta_tau;
u_min_final = u_min_lower_limit+y*delta_u_min;

disp(tau_final);
disp(u_min_final);

% Defining the spectrum of u_min's values
u_min_lower_limit = u_min * 0.5;
u_min_upper_limit = u_min * 1.5;
u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);

% Defining the spectrum of tau's values
tau_lower_limit = tau * 0.5;
tau_upper_limit = tau * 1.5;
tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);
delta_chi = chi_mat - chi_mat(y, x);
[u_min_axis, tau_axis] = meshgrid(u_min_vec, tau_vec);
contour(u_min_axis, tau_axis, delta_chi, [1 2.3 2.71 6.63],'ShowText', true)

% Set axes origin at tau_final & u_min_final.
% Calculate the projaction from the origin to the furtherest point on the
% elipse, showing different contour values of chi-squared at points 1 2.3 2.71 6.63

