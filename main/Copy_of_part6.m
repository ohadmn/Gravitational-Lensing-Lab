function [tau_final, u_min_final, fbl_final, T0_final] = Copy_of_part6(file_path, T0, fbl, tau, u_min)
    % Load and process the data
    data = load(file_path);
    time = data(:,1);
    magnitude = data(:,2);
    err = data(:,3);
    
    % Convert magnitude to flux
    flux = 10.^(0.4 * (-magnitude));
    sigma_flux = 0.4 * log(10) * err .* flux;
    
    % Normalize flux
    f_norm = flux(1);
    flux = flux / f_norm;
    sigma_flux = sigma_flux / f_norm;       

    % Define the chi-squared function
    function chi2 = chi_squared(params)
        tau = params(1);
        u_min = params(2);
        fbl = params(3);
        T0 = params(4);

        delta_t = time - T0;
        t = delta_t / tau;
        u = sqrt(u_min^2 + t.^2);
        mu = (2 + u.^2) ./ (u .* sqrt(4 + u.^2));
        new_mu = (1 - fbl) + fbl .* mu;
        chi2 = sum(((flux - new_mu) ./ sigma_flux).^2);
    end

    % Initial guesses for the parameters
    T0_initial = T0 * 0.99;        
    fbl_initial = fbl * 0.8;              
    tau_initial = tau * 0.8;
    u_min_initial = u_min * 0.8;  

    % Initial parameter guesses
    initial_params = [tau_initial, u_min_initial, fbl_initial, T0_initial];

    % Minimize the chi-squared function
    [final_params] = fminsearch(@chi_squared, initial_params);

    % Extract the optimized parameters
    tau_final = final_params(1);
    u_min_final = final_params(2);
    fbl_final = final_params(3);
    T0_final = final_params(4);
end
