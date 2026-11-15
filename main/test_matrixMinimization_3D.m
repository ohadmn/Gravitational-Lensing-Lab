numIterations = 10000; 
dim = 500; 
relativeDeviations = zeros(numIterations, 1); % Store relative deviations

for iter = progress(1:numIterations)
    [X, Y, Z] = ndgrid(1:dim, 1:dim, 1:dim);
    mat = abs(dim/2 - X) + abs(dim/2 - Y) + abs(dim/2 - Z) + randi([1, max(X(:)) + max(Y(:)) + max(Z(:))], size(X));

    % Find minimum using our algorithm
    [i_minimal, j_minimal, k_minimal] = matrixMinimization_3D(mat);
    myMinValue = mat(i_minimal, j_minimal, k_minimal);

    % Find minimum using MATLAB's built-in function for comparison
    [minValue] = min(mat(:));

    % Calculate the range of matrix values
    rangeMat = max(mat(:)) - minValue;

    % Calculate relative deviation as a percentage
    relativeDeviations(iter) = abs(myMinValue - minValue) / rangeMat * 100;
end

meanRelativeDeviation = mean(relativeDeviations);
maxRelativeDeviation = max(relativeDeviations);
minRelativeDeviation = min(relativeDeviations);
stdRelativeDeviation = std(relativeDeviations);

disp('Statistical Summary of Relative Deviation (as percentage):');
disp(['Mean Relative Deviation: ', num2str(meanRelativeDeviation), '%']);
disp(['Max Relative Deviation: ', num2str(maxRelativeDeviation), '%']);
disp(['Min Relative Deviation: ', num2str(minRelativeDeviation), '%']);
disp(['Standard Deviation of Relative Deviation: ', num2str(stdRelativeDeviation), '%']);
