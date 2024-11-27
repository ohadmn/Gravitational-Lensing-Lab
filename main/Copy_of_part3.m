function[tau_final, u_min_final, T0_final] = Copy_of_part3(file_path, T0, tau, u_min)
    % Hardcoded data per incident.
    data = load(file_path);
    time = data(:,1);
    magnitude = data(:,2);
    % Converting magnitude to flux based on eq.(1.14) when F0=1 and m0=0
    flux = 1*10.^(0.4*(-magnitude));
    err = data(:,3);
    % Converting magnitude to flux based on the derivative of eq.(1.14) with respect to m, when F0=1 and m0=0
    sigma_flux = -0.4*log(10).*err.*10.^(0.4*(-magnitude));
    % T0 = T0;
    
    % Normalizing the flux values with respect to flux value without
    % microlensing.
    f_norm = flux(1);
    flux = flux/f_norm;
    sigma_flux = sigma_flux/(f_norm);
    
    num_measurements = 20;
    
    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min*0.5;
    u_max=u_min*1.5;
    delta_u_min=(u_max-u_min_lower_limit)/num_measurements;
    %delta_u_min = (u_min*0.5)/(num_measurements/2);
    
    % Defining the spectrum of tau's values
    tau_lower_limit = tau*0.5;
    tau_max=tau*1.5;
    %delta_tau = (tau*0.5)/(num_measurements/2);
    delta_tau=(tau_max-tau_lower_limit)/num_measurements;
    
    % Defining the spectrum of T0's values
    T0_lower_limit = T0*0.999999;
    T0_max = T0*1.000001;
    % 
    % T0_lower_limit = T0-10;
    % T0_max=T0+10;
    % delta_T0 = (T0*0.5)/(num_measurements/2);
    delta_T0=(T0_max-T0_lower_limit)/num_measurements;
    
    
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements, num_measurements, num_measurements);
    flux_size = length(flux);
    % Calculating chi squared with different variations from tau and u_min
    % spectrum.
    for l = 1:num_measurements
        for i = 1:num_measurements
            for j = 1:num_measurements           
                for k = 1:flux_size
                    t = (time(k)-(T0_lower_limit + delta_T0 * l))/(tau_lower_limit + delta_tau * i);
                    u = ((u_min_lower_limit + delta_u_min * j) ^ 2 + t ^ 2) ^ (1/2);
                    mu = (2 + u ^ 2) / (u * (4 + u ^ 2) ^ (1/2));
                    
                    chi_mat(i, j, l) = chi_mat(i, j, l) + (((flux(k) - mu) / sigma_flux(k)) ^ 2);
                end
            end
        end
    end
    % Finding the minimazation of the matrix.
    %[x, y] = matrixMinimization(chi_mat);
    % Test - Compare with MATLAB built-in functions
    [minValue, linearIndex] = min(chi_mat(:)); % Find the minimum value and its linear index
    [y, x, z] = ind2sub(size(chi_mat), linearIndex); % Convert linear index to row and column indices
    disp(x);
    disp(y);
    disp(z);
    
    
    tau_final = tau_lower_limit + x * delta_tau;
    u_min_final = u_min_lower_limit + y * delta_u_min;
    T0_final = T0_lower_limit + z * delta_T0;

    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min * 0.5;
    u_min_upper_limit = u_min * 1.5;
    u_min_vec = linspace(u_min_lower_limit, u_min_upper_limit, num_measurements);
    
    % Defining the spectrum of tau's values
    tau_lower_limit = tau * 0.5;
    tau_upper_limit = tau * 1.5;
    tau_vec = linspace(tau_lower_limit, tau_upper_limit, num_measurements);

    % Defining the spectrum of tau's values
    T0_lower_limit = T0 * 0.8;
    T0_upper_limit = T0 * 1.2;
    T0_vec = linspace(T0_lower_limit, T0_upper_limit, num_measurements);    

    delta_chi = chi_mat - chi_mat(y, x, z);
    delta_chi(3) = delta_chi(3 ).*10^-6 ; 

    delta_chi_tau_u_min = delta_chi(:, :, z);
    figure;
    % subplot(1, 3, 1);
    contour(u_min_vec, tau_vec, delta_chi_tau_u_min, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('u_{min} (unitless)')
    ylabel('\tau')
    title('Contour of u_{min} - \tau plane')
    
    delta_chi_T_0_u_min = squeeze(delta_chi(y, :, :));
    figure;
    % subplot(1, 3, 2);
    contour(T0_vec, u_min_vec, delta_chi_T_0_u_min, [1 2.3 2.71 6.63] ,'ShowText', true)
    xlabel('\T_0 (unitless)')
    ylabel('u_{min}')
    title('Contour of u_min - T_0 plane')
    
    delta_chi_tau_T_0 = squeeze(delta_chi(:, x, :));
    figure;
    % subplot(1, 3, 3);
    contour(T0_vec, tau_vec, delta_chi_tau_T_0, [1 2.3 2.71 6.63],'ShowText', true)
    xlabel('\T_0 (unitless)')
    ylabel('\tau')
    title('Contour of \T_0 - \tau plane')
    

end