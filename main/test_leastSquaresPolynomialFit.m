% Define true values of the polynomial coefficients for the test
a_0 = 1;
b_0 = 2;
c_0 = 3;

% Initialize arrays to store computed values of a, b, c, and chi for each iteration
a_val = [];
b_val = [];
c_val = [];
chi_val = [];

% Run the loop 10,000 times to test the robustness of the fitting function
for i = 1:10000
    % Generate 100 equally spaced values of x between 1 and 5
    x = linspace(1,5,100);
    
    % Calculate the true y values without noise
    y_org = (a_0 .* x.^2 + b_0 .* x + c_0);
    
    % Add random Gaussian noise to the y values
    y = y_org + 0.01 * randn(1,100);
    
    % Define a constant sigma array (all values are 0.1) representing the uncertainties
    sigma = 0.1 * ones(1,100);
    
    % Estimate coefficients a, b, and c using the fitting function parab
    [a, b, c] = parab(x, y, sigma);
    
    % Store the computed values of a, b, and c for this iteration
    a_val(i) = a;
    b_val(i) = b;
    c_val(i) = c;
    
    % Calculate chi-squared (goodness of fit) and store it
    chi_val(i) = sum((y - a .* x.^2 - b .* x - c).^2 ./ sigma.^2);
end

% Calculate the mean value of each coefficient and its standard deviation.
mean_a = mean(a_val);
std_a = std(a_val);

mean_b = mean(b_val);
std_b = std(b_val);

mean_c = mean(c_val);
std_c = std(c_val);

mean_chi = mean(chi_val);
std_chi = std(chi_val);

fprintf('Estimated value of a: %.3f ± %.3f\n', mean_a, std_a);
fprintf('Estimated value of b: %.3f ± %.3f\n', mean_b, std_b);
fprintf('Estimated value of c: %.3f ± %.3f\n', mean_c, std_c);
fprintf('Estimated value of c: %.3f ± %.3f\n', mean_chi, std_chi);


% Plot histograms of the computed values of a, b, and c, and chi-squared

% Histogram for coefficient a
subplot(2, 2, 1);
histogram(a_val)
xlabel('Iterations (unitless)')
ylabel('Count')
title('Histogram of a (Estimated Coefficient)')

% Histogram for coefficient b
subplot(2, 2, 2);
histogram(b_val)
xlabel('Iterations (unitless)')
ylabel('Count')
title('Histogram of b (Estimated Coefficient)')

% Histogram for coefficient c
subplot(2, 2, 3);
histogram(c_val)
xlabel('Iterations (unitless)')
ylabel('Count')
title('Histogram of c (Estimated Coefficients)')

% Histogram for chi-squared values
subplot(2, 2, 4);
histogram(chi_val)
xlabel('Iterations (unitless)')
ylabel('Count')
title('Histogram of \chi^2 (Goodness of Fit)')
