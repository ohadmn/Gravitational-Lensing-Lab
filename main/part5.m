function[tau_final, u_min_final, fbl_final] = part5(file_path, T0, fbl, tau, u_min)
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
    
    num_measurements = 10;
    
    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min*0.5;
    delta_u_min = (u_min*0.5)/(num_measurements/2);
    % Defining the spectrum of tau's values
    tau_lower_limit = tau*0.5;
    delta_tau = (tau*0.5)/(num_measurements/2);
    % Defining the spectrum of T0's values
    T0_lower_limit = T0*0.5;
    delta_T0 = (T0*0.5)/(num_measurements/2);
    % Defining the spectrum of fbl's values
    fbl_lower_limit = fbl*0.5;
    delta_fbl = (fbl*0.5)/(num_measurements/2);
    
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements, num_measurements, num_measurements, num_measurements);
    flux_size = length(flux);
    % Calculating chi squared with different variations from tau and u_min
    % spectrum.
    for l = 1:num_measurements
        for i = 1:num_measurements
            for j = 1:num_measurements  
                for n = 1:num_measurements
                    for k = 1:flux_size
                        t = (time(k)-(T0_lower_limit + delta_T0 * n))/(tau_lower_limit + delta_tau * i);
                        u = ((u_min_lower_limit + delta_u_min * j) ^ 2 + t ^ 2) ^ (1/2);
                        mu = (2 + u ^ 2) / (u * (4 + u ^ 2) ^ (1/2));
                        fbl = fbl_lower_limit + delta_fbl * l;
                        new_mu = 1 - fbl + (fbl * mu);
                        chi_mat(i, j, l, n) = chi_mat(i, j, l, n) + (((flux(k) - new_mu) / sigma_flux(k)) ^ 2);
                    end
                end
            end
        end
    end
    % Finding the minimazation of the matrix.
    %[x, y] = matrixMinimization(chi_mat);
    % Test - Compare with MATLAB built-in functions
    [minValue, linearIndex] = min(chi_mat(:)); % Find the minimum value and its linear index
    [i, j, l, n] = ind2sub(size(chi_mat), linearIndex); % Convert linear index to row and column indices
    disp(j);
    disp(i);
    disp(l);
    disp(n);
    
    tau_final = tau_lower_limit + i * delta_tau;
    u_min_final = u_min_lower_limit + j * delta_u_min;
    fbl_final = fbl_lower_limit + l * delta_fbl;
    T0_final = T0_lower_limit + n * delta_T0;


    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min * 0.5;
    u_min_upper_limit = u_min * 1.5;
    u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);
    
    % Defining the spectrum of tau's values
    tau_lower_limit = tau * 0.5;
    tau_upper_limit = tau * 1.5;
    tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);

    % Defining the spectrum of T0's values
    T0_lower_limit = T0 * 0.8;
    T0_upper_limit = T0 * 1.2;
    T0_vec = linspace(T0_lower_limit, T0_upper_limit, num_measurements);
    
    % Defining the spectrum of fbl's values
    fbl_lower_limit = fbl * 0.8;
    fbl_upper_limit = fbl * 1.2;
    fbl_vec = linspace(fbl_lower_limit, fbl_upper_limit, num_measurements);

    delta_chi = chi_mat - chi_mat(j, i, l, n);
    
    delta_chi_tau_u_min = squeeze(delta_chi(:, :, l, n));
    figure;
    % subplot(1, 6, 1);
    contour(u_min_vec, tau_vec, delta_chi_tau_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('\tau')
    title('Contour of u_{min} - \tau plane')
    
    delta_chi_fbl_u_min = squeeze(delta_chi(i, :, :, n));
    figure;
    % subplot(1, 6, 2);
    contour(fbl_vec, u_min_vec, delta_chi_fbl_u_min, [1 2.3 2.71 6.63] ,'ShowText', true)
    xlabel('fbl (unitless)')
    ylabel('u_{min}')
    title('Contour of u_{min} - fbl plane')

    delta_chi_T0_u_min = squeeze(delta_chi(j, :, l, :));
    figure;
    % subplot(1, 6, 3);
    contour(T0_vec, u_min_vec, delta_chi_T0_u_min, [1 2.3 2.71 6.63] ,'ShowText', true)
    xlabel('T_0 (unitless)')
    ylabel('u_{min}')
    title('Contour of T_0 - u_{min} plane')

    delta_chi_tau_fbl = squeeze(delta_chi(:, i, :, n));
    figure;
    % subplot(1, 6, 4);
    contour(fbl_vec, tau_vec, delta_chi_tau_fbl, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('fbl (unitless)')
    ylabel('\tau')
    title('Contour of fbl - \tau plane')

    delta_chi_tau_T0 = squeeze(delta_chi(:, i, l, :));
    figure;
    % subplot(1, 6, 5);
    contour(T0_vec, tau_vec, delta_chi_tau_T0, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('T_0 (unitless)')
    ylabel('\tau')
    title('Contour of T_0 - \tau plane')

    delta_chi_fbl_T0 = squeeze(delta_chi(j, i, :, :));
    figure;
    % subplot(1, 6, 6);
    contour(T0_vec, fbl_vec, delta_chi_fbl_T0, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('T_0 (unitless)')
    ylabel('fbl')
    title('Contour of T_0 - fbl plane')
    

    % The distance of the contour elipses from the origin is the delta of the answer.
    % Calculate the projection of each contour from the origin.
    % The origin is defined by the vaules returned from this function(minimum).
    % tau_final\ fbl_final \ u_min_final 
    

end