function[tau_final, u_min_final, fbl_final, T0_final] = Copy_of_part5(file_path, T0, fbl, tau, u_min)
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
    
    num_measurements = 1000;
    
   % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min * 0.5;
    u_min_upper_limit = u_min * 1.5;
    delta_u_min = (u_min_upper_limit - u_min_lower_limit) / num_measurements;
    % Defining the spectrum of tau's values
    tau_lower_limit = tau * 0.5;
    tau_upper_limit = tau * 1.5;
    delta_tau = (tau_upper_limit - tau_lower_limit) / num_measurements;
    % Defining the spectrum of T0's values
    T0_lower_limit = T0 * 0.999999;
    T0_upper_limit = T0 * 1.000001;
    delta_T0 = (T0_upper_limit - T0_lower_limit) / num_measurements;
    % Defining the spectrum of fbl's values
    fbl_lower_limit = fbl * 0.9;
    fbl_upper_limit = fbl * 1.1;
    delta_fbl = (fbl_upper_limit - fbl_lower_limit) / num_measurements;
    
    % Defining the spectrum of u_min's values
    u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);
    
    % Defining the spectrum of tau's values
    tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);

    % Defining the spectrum of T0's values
    T0_vec = linspace(T0_lower_limit, T0_upper_limit, num_measurements);
    
    % Defining the spectrum of fbl's values
    fbl_vec = linspace(fbl_lower_limit, fbl_upper_limit, num_measurements);
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements, num_measurements);
     
    % 2d matrix of u_min and tau
    % Precompute delta_time
    delta_time = time - T0;  % size (574,1)
    % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        % Compute t for all tau_values
        t = delta_time(k) ./ tau_vec;  % size (1, num_measurements)        
        % Create meshgrid for u_min_values and t
        [u_min_grid, t_grid] = meshgrid(u_min_vec, t);
        % Compute u
        u = sqrt(u_min_grid.^2 + t_grid.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl + (fbl * mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end

    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_tau_u_min = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(u_min_vec, tau_vec, delta_chi_tau_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('\tau')
    title('Contour of u_{min} - \tau plane')
    
    % 2d matrix of u_min and fbl
    % Precompute delta_time
    delta_time = time - T0;  % size (574,1)
    % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        % Compute t for all tau_values
        t = delta_time(k) / tau;  % size (1, num_measurements)
        [u_min_grid, fbl_grid] = meshgrid(u_min_vec, fbl_vec);
        % Compute u
        u = sqrt(u_min_grid.^2 + t.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl_grid + (fbl_grid .* mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end

    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_fbl_u_min = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(u_min_vec, fbl_vec, delta_chi_fbl_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('fbl')
    title('Contour of u_{min} - fbl plane')

    % 2d matrix of u_min and T0
    % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        [u_min_grid, T0_grid] = meshgrid(u_min_vec, T0_vec);
        delta_time = time(k) - T0_grid;  % size (574,1)
        % Compute t for all tau_values
        t = delta_time / tau;  % size (1, num_measurements)
        % Compute u
        u = sqrt(u_min_grid.^2 + t.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl + (fbl * mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end

    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_T0_u_min = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(u_min_vec, T0_vec, delta_chi_T0_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('T0')
    title('Contour of u_{min} - T0 plane')
    
    % 2d matrix of T0 and tau
    % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        [tau_grid, T0_grid] = meshgrid(tau_vec, T0_vec);
        delta_time = time(k) - T0_grid;  % size (574,1)
        % Compute t for all tau_values
        t = delta_time ./ tau_grid;  % size (1, num_measurements)
        % Compute u
        u = sqrt(u_min^2 + t.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl + (fbl * mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end

    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_T0_tau = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(tau_vec, T0_vec, delta_chi_T0_tau, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('\tau (unitless)')
    ylabel('T0')
    title('Contour of \tau - T0 plane')

    % % 2d matrix of fbl and tau
    % % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        [tau_grid, fbl_grid] = meshgrid(tau_vec, fbl_vec);
        delta_time = time(k) - T0;  % size (574,1)
        % Compute t for all tau_values
        t = delta_time ./ tau_grid;  % size (1, num_measurements)
        % Compute u
        u = sqrt(u_min^2 + t.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl_grid + (fbl_grid .* mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end

    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_fbl_tau = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(tau_vec, fbl_vec, delta_chi_fbl_tau, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('\tau (unitless)')
    ylabel('fbl')
    title('Contour of \tau - fbl plane')

    % 2d matrix of fbl and T0
    % Loop over k to accumulate chi-squared values
    for k = progress(1:length(time))
        [T0_grid, fbl_grid] = meshgrid(T0_vec, fbl_vec);
        delta_time = time(k) - T0_grid;  % size (574,1)
        % Compute t for all tau_values
        t = delta_time / tau;  % size (1, num_measurements)
        % Compute u
        u = sqrt(u_min^2 + t.^2);  % size (num_measurements, num_measurements)
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = 1 - fbl_grid + (fbl_grid .* mu);
        % Compute residuals
        residual = ((flux(k) - new_mu) / sigma_flux(k)).^2;
        % Accumulate chi-squared values
        chi_mat = chi_mat + residual;
    end
    
    [x, y] = matrixMinimization_2D(chi_mat);
    delta_chi_fbl_T0 = chi_mat - chi_mat(x, y);
    figure;
    % subplot(1, 3, 1);
    contour(T0_vec, fbl_vec, delta_chi_fbl_T0, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('T0 (unitless)')
    ylabel('fbl')
    title('Contour of T0 - fbl plane')
end