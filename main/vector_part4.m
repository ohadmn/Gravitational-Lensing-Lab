function[tau_final, u_min_final, fbl_final] = vector_part4(file_path, T0, fbl, tau, u_min)
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
    sigma_flux = sigma_flux / (f_norm);
    
    num_measurements = 600;
    
    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min * 0.5;
    u_min_upper_limit = u_min * 1.5;
    delta_u_min = (u_min_upper_limit - u_min_lower_limit) / num_measurements;
    u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);
    % Defining the spectrum of tau's values
    tau_lower_limit = tau * 0.5;
    tau_upper_limit = tau * 1.5;
    delta_tau = (tau_upper_limit - tau_lower_limit) / num_measurements;
    tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);
    % Defining the spectrum of fbl's values
    fbl_lower_limit = fbl * 0.88;
    fbl_upper_limit = fbl * 1.12;
    delta_fbl = (fbl_upper_limit - fbl_lower_limit) / num_measurements;
    fbl_vec = linspace(fbl_lower_limit, fbl_upper_limit, num_measurements);
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements, num_measurements, num_measurements);
    [tau_grid, u_min_grid, fbl_grid] = ndgrid(tau_vec, u_min_vec, fbl_vec);

    % Calculating chi squared with different variations from tau and u_min
    % spectrum.
    for k = progress(1:length(time))
        delta_t = time(k) - T0;
        % Compute t for all combinations of tau
        t = delta_t ./ tau_grid;
        
        % Compute u for all combinations of u_min and t
        u = sqrt(u_min_grid.^2 + t.^2);
        
        % Compute mu
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        
        % Compute new_mu for all combinations of mu and fbl
        new_mu = 1 - fbl_grid + fbl_grid .* mu;
        
        % Compute chi-squared values
        chi_mat = chi_mat + ((flux(k) - new_mu) / sigma_flux(k)).^2;
    end

    % Finding the minimazation of the matrix.
    [x, y, z] = matrixMinimization_3D(chi_mat);
    disp(x);
    disp(y);
    disp(z);

    tau_final = tau_lower_limit + x * delta_tau;
    u_min_final = u_min_lower_limit + y * delta_u_min;
    fbl_final = fbl_lower_limit + z * delta_fbl;

    delta_chi = chi_mat - chi_mat(y, x, z);

    delta_chi_tau_u_min = delta_chi(:, :, z);
    [u_min_axis, tau_axis] = meshgrid(u_min_vec, tau_vec);    
    figure;
    % subplot(1, 3, 1);
    contour(u_min_axis, tau_axis, delta_chi_tau_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('\tau')
    title('Contour of u_{min} - \tau plane')
    
    delta_chi_fbl_u_min = squeeze(delta_chi(y, :, :));
    [fbl_axis, u_min_axis] = meshgrid(fbl_vec, u_min_vec);
    figure;
    % subplot(1, 3, 2);
    contour(fbl_axis, u_min_axis, delta_chi_fbl_u_min, [1 2.3 2.71 6.63] ,'ShowText', true)
    xlabel('fbl (unitless)')
    ylabel('u_{min}')
    title('Contour of u_{min} - fbl plane')
    
    delta_chi_tau_fbl = squeeze(delta_chi(:, x, :));
    [fbl_axis, tau_axis] = meshgrid(fbl_vec, tau_vec);
    figure;
    % subplot(1, 3, 3);
    contour(fbl_axis, tau_axis, delta_chi_tau_fbl, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('fbl (unitless)')
    ylabel('\tau')
    title('Contour of fbl - \tau plane')

    % The distance of the contour elipses from the origin is the delta of the answer.
    % Calculate the projection of each contour from the origin.
    % The origin is defined by the vaules returned from this function(minimum).
    % tau_final\ fbl_final \ u_min_final 
    

end