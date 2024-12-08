function [x, y] = matrixMinimization_2D(mat)
    % Initialize variables
    [dimX, dimY] = size(mat);
    currentMin = flintmax; % Initialize to the largest possible value
    x = 1; y = 1;
    minUpdated = true;

    % Loop until no smaller value is found
    while minUpdated
        minUpdated = false;

        % Search along the x-dimension (fix y)
        for i = 1:dimX
            if mat(i, y) < currentMin
                currentMin = mat(i, y);
                x = i;
                minUpdated = true;
            end
        end

        % Search along the y-dimension (fix x)
        for j = 1:dimY
            if mat(x, j) < currentMin
                currentMin = mat(x, j);
                y = j;
                minUpdated = true;
            end
        end
    end

    % Return the indices of the minimum value
end
