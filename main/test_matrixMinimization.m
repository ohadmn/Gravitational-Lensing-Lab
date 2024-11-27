mat = ones(1000,1000);
mat = mat + 0.1*randn(1000);

[u_min_index, tau_index] = matrixMinimization(mat);

% Test - Compare with MATLAB built-in functions
[minValue, linearIndex] = min(mat(:)); % Find the minimum value and its linear index
[row, col] = ind2sub(size(mat), linearIndex); % Convert linear index to row and column indices

% Display results
disp(['Minimum Value: ', num2str(minValue)]);
disp(['Row Index: ', num2str(row)]);
disp(['Column Index: ', num2str(col)]);