function [x, y] = matrixMinimization_2D(mat)
    % Initialize variables
    [rows, cols] = size(mat);
    currentMin = flintmax;
    row = 1;
    col = 1;
    minUpdated = true;
    searchDirection = 1; % 1 for row search, 2 for column search

    % Loop until no smaller value is found
    while minUpdated
        minUpdated = false;

        if searchDirection == 1
            % Search along the current row
            for c = 1:cols
                if mat(row, c) < currentMin
                    currentMin = mat(row, c);
                    col = c;
                    minUpdated = true;
                end
            end
            searchDirection = 2; % Switch to column search
        else
            % Search along the current column
            for r = 1:rows
                if mat(r, col) < currentMin
                    currentMin = mat(r, col);
                    row = r;
                    minUpdated = true;
                end
            end
            searchDirection = 1; % Switch to row search
        end
    end

    % Return the coordinates of the global minimum
    x = col;
    y = row;
end
