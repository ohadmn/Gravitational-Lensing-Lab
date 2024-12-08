dim = 1000;
mat  = zeros(dim, dim);
for i = 1:dim
    for j = 1:dim
        mat(i, j) = abs(dim/2 - i) + abs(dim/2 - j) + abs(randi([1 i+j], 1));
    end
end

% Our minimization results
[i_minimal, j_minimal] = matrixMinimization_2D(mat);
disp("matrixMinimization_2D:")
disp(['Minimum Value: ', num2str(mat(i_minimal, j_minimal))]);
disp(['Row Index: ', num2str(i_minimal)]);
disp(['Column Index: ', num2str(j_minimal)]);

% Test - Compare with MATLAB built-in functions
[minValue, linearIndex] = min(mat(:)); % Find the minimum value and its linear index
[row, col] = ind2sub(size(mat), linearIndex); % Convert linear index to row and column indices

% Display Matlab's results
disp("Matlab's built-in function:")
disp(['Minimum Value: ', num2str(minValue)]);
disp(['Row Index: ', num2str(row)]);
disp(['Column Index: ', num2str(col)]);