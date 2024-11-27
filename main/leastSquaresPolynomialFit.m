% leastSquaresPolynomialFit.m
function [max_time, flux_max] = leastSquaresPolynomialFit(filePath)
    % Load the data file
    data = load(filePath);
    time = data(:,1);
    magnitude = data(:,2);
    err = data(:,3);

    % Convert magnitude to flux
    flux = 1 * 10.^(0.4 * (-magnitude));
    sigma_flux = -0.4 * log(10) * err .* 10.^(0.4 * (-magnitude));

    % Normalize the flux values
    flux_norm_factor = (5 * 10^-9);
    flux = flux / flux_norm_factor;
    sigma_flux = sigma_flux / flux_norm_factor;

    % Slicing interval from OGLE's data
    flux_sliced_values = flux(493:533);         % Extract 41 values from flux
    time_sliced_values = time(493:533);         % Extract corresponding time values
    sigma_sliced_values = sigma_flux(493:533);  % Extract corresponding error values

    % Normalize time with respect to the first value
    time_norm_factor = time_sliced_values(1);
    time_sliced_values = time_sliced_values - time_norm_factor;

    % Calculate a, b, c using the parab function
    [a, b, c] = parab(time_sliced_values, flux_sliced_values, sigma_sliced_values);

    % Define the parabola fit
    parabola = a * time_sliced_values.^2 + b * time_sliced_values + c;

    % Calculate the time and flux of the maximum value
    max_time = -b / (2 * a) + time_norm_factor;
    flux_max = a * (max_time - time_norm_factor)^2 + b * (max_time - time_norm_factor) + c;

    % Display results
    disp("Time of max flux (Julian day units):");
    disp(max_time);
    disp("Max flux:");
    disp(flux_max);

    % Plot the results
    figure;
    plot(time_sliced_values, parabola, '-');
    hold on;
    errorbar(time_sliced_values, flux_sliced_values, sigma_sliced_values, ".", 'MarkerSize', 10);
    hold off;
    title('Flux as a function of time (i = 493:533)');
    xlabel('Time (T_{0})');
    ylabel('Flux (5 * 10^-9)');
end
