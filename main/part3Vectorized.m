function[tau_final, u_min_final, T0_final] = part3Vectorized(file_path, T0, tau, u_min)
    % Hardcoded data per incident.
    data = load(file_path);
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
    
    num_measurements = 600;
    
    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min * 0.5;
    u_min_upper_limit = u_min * 1.5;
    delta_u_min = (u_min_upper_limit - u_min_lower_limit) / num_measurements;
    
    % Defining the spectrum of tau's values
    tau_lower_limit = tau * 0.5;
    tau_upper_limit = tau * 1.5;
    delta_tau=(tau_upper_limit - tau_lower_limit) / num_measurements;
    
    % Defining the spectrum of T0's values
    T0_lower_limit = T0 * 0.999999;
    T0_upper_limit = T0 * 1.000001;
    delta_T0 = (T0_upper_limit - T0_lower_limit) / num_measurements;

    % Defining the spectrum of u_min's values
    u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);
    
    % Defining the spectrum of tau's values
    tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);

    % Defining the spectrum of tau's values
    T0_vec = linspace(T0_lower_limit, T0_upper_limit, num_measurements);
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements, num_measurements, num_measurements);
    % Calculating chi squared with different variations from tau and u_min
    % spectrum.
    
    % Create parameter grids using ndgrid
    [tau_grid, u_min_grid, T0_grid] = ndgrid(tau_vec, u_min_vec, T0_vec);

    % Loop over data points
    for k = progress(1:length(time))
        % Compute t for all combinations of tau and T0
        delta_t = time(k) - T0_grid;  % delta_t is (num_measurements x num_measurements x num_measurements)
        t = delta_t ./ tau_grid;

        % Compute u
        u = sqrt(u_min_grid.^2 + t.^2);

        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));

        % Compute residuals and accumulate chi-squared values
        residual = ((flux(k) - mu) / sigma_flux(k)).^2;
        chi_mat = chi_mat + residual;
    end

    % Finding the minimazation of the matrix.
    [x, y, z] = matrixMinimization_3D(chi_mat);
    disp(x);
    disp(y);
    disp(z);
    
    tau_final = tau_lower_limit + x * delta_tau;
    u_min_final = u_min_lower_limit + y * delta_u_min;
    T0_final = T0_lower_limit + z * delta_T0;    

    delta_chi = chi_mat - chi_mat(y, x, z);

    delta_chi_tau_u_min = delta_chi(:, :, z);
    figure;
    subplot(1, 3, 1);
    contour(u_min_vec, tau_vec, delta_chi_tau_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('\tau')
    title('Contour of u_{min} - \tau plane')
    
    delta_chi_T_0_u_min = squeeze(delta_chi(y, :, :));
    %figure;
    subplot(1, 3, 2);
    contour(T0_vec, u_min_vec, delta_chi_T_0_u_min, [1 2.3 2.71 6.63] ,'ShowText', true)
    xlabel('T_0 (unitless)')
    ylabel('u_{min}')
    title('Contour of u_{min} - T_0 plane')
    
    delta_chi_tau_T_0 = squeeze(delta_chi(:, x, :));
    %figure;
    subplot(1, 3, 3);
    contour(T0_vec, tau_vec, delta_chi_tau_T_0, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('T_0 (unitless)')
    ylabel('\tau')
    title('Contour of T_0 - \tau plane')
    

end