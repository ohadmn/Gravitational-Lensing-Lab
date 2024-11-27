function[tau_final] = part2(file_path, T0, tau, u_min)
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
    
    num_measurements = 100;
    
    % Defining the spectrum of u_min's values
    u_min_lower_limit = u_min*0.5;
    delta_u_min = (u_min*0.5)/(num_measurements/2);
    % Defining the spectrum of tau's values
    tau_lower_limit = tau*0.5;
    delta_tau = (tau*0.5)/(num_measurements/2);
    
    % The dimensions of chi_mat are the number of measurements.
    chi_mat = zeros(num_measurements,num_measurements);
    flux_size = length(flux);
    % Calculating chi squared with different variations from tau and u_min
    % spectrum.
    for i = 1:num_measurements
        for j = 1:num_measurements           
            for k = 1:flux_size
                t = (time(k)-T0)/(tau_lower_limit + delta_tau * i);
                u = ((u_min_lower_limit + delta_u_min * j) ^ 2 + t ^ 2) ^ (1/2);
                mu = (2 + u ^ 2) / (u * (4 + u ^ 2) ^ (1/2));
                
                chi_mat(i,j) = chi_mat(i,j) + (((flux(k) - mu) / sigma_flux(k)) ^ 2);
            end
        end
    end
    
    % Finding the minimazation of the matrix.
    [x, y] = matrixMinimization(chi_mat);
    % % Test - Compare with MATLAB built-in functions
    % [minValue, linearIndex] = min(chi_mat(:)); % Find the minimum value and its linear index
    % [y, x] = ind2sub(size(chi_mat), linearIndex); % Convert linear index to row and column indices
    disp(x);
    disp(y);
    
    
    tau_final = tau_lower_limit+x*delta_tau;
    u_min_final = u_min_lower_limit+y*delta_u_min;
    

end