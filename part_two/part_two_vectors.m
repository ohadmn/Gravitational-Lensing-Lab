% Hardcoded data per incident.
data = load("~/Desktop/phot.dat");
time = data(:,1);
magnitude = data(:,2);
% Converting magnitude to flux based on eq.(1.14) when F0=1 and m0=0
flux = 1*10.^(0.4*(-magnitude));
err = data(:,3);
% Converting magnitude to flux based on the derivative of eq.(1.14) with respect to m, when F0=1 and m0=0
sigma_flux = -0.4*log(10).*err.*10.^(0.4*(-magnitude));

% Normalizing the flux values with respect to flux value without microlensing.
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
u_min_values = u_min_lower_limit + delta_u_min * (1:num_measurements);

% Defining the spectrum of tau's values
tau_lower_limit = tau*0.5;
delta_tau = (tau*0.5)/(num_measurements/2);
tau_values = tau_lower_limit + delta_tau * (1:num_measurements);

% Pre-allocate chi_mat
chi_mat = zeros(num_measurements, num_measurements);

% Precompute delta_time
delta_time = time - T0;  % size (574,1)

% Loop over k to accumulate chi-squared values
for k = 1:length(time)
    % Compute t for all tau_values
    t = delta_time(k) ./ tau_values;  % size (1, num_measurements)

    % Create meshgrid for u_min_values and t
    [u_min_grid, t_grid] = meshgrid(u_min_values, t);

    % Compute u
    u = sqrt(u_min_grid.^2 + t_grid.^2);  % size (num_measurements, num_measurements)

    % Compute mu
    mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));

    % Compute residuals
    residual = ((flux(k) - mu) / sigma_flux(k)).^2;

    % Accumulate chi-squared values
    chi_mat = chi_mat + residual;
end

% Finding the minimization of the matrix.
[x, y] = matrix_minimization(chi_mat);
disp(x);
disp(y);

tau_final = tau_lower_limit + x * delta_tau;
u_min_final = u_min_lower_limit + y * delta_u_min;

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
