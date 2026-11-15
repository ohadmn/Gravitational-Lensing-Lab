function [minCol, minRow] = matrix_mini_debug(matrix)
    % Ensure the matrix is square
    assert(size(matrix, 1) == size(matrix, 2), 'Matrix must be square!');
    
    % Initialize variables
    [matrixSize, ~] = size(matrix); % Dimensions of the square matrix
    currentMin = flintmax; % Start with the largest possible value
    currentRow = 1;
    currentCol = 1;

    % Iterate through the entire matrix
    for row = 1:matrixSize
        for col = 1:matrixSize
            % Check if the current element is smaller than the current minimum
            if matrix(row, col) < currentMin
                currentMin = matrix(row, col);
                currentRow = row;
                currentCol = col;
            end
        end
    end

    % Output the indices of the global minimum
    minRow = currentRow;
    minCol = currentCol;

    % Debugging print
    fprintf('Row = %d, Column = %d, Minimum Value = %.5f\n', ...
            minRow, minCol, currentMin);
end
