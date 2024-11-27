function [x, y, z] = matrixMinimization_3D(mat)
    % Initialize variables
    [dimX, dimY, dimZ] = size(mat);
    currentMin = flintmax; % Initialize to the largest possible value
    x = 1; y = 1; z = 1;
    minUpdated = true;

    % Loop until no smaller value is found
    while minUpdated
        minUpdated = false;

        % Search along the x-dimension (fix y and z)
        for i = 1:dimX
            if mat(i, y, z) < currentMin
                currentMin = mat(i, y, z);
                x = i;
                minUpdated = true;
            end
        end

        % Search along the y-dimension (fix x and z)
        for j = 1:dimY
            if mat(x, j, z) < currentMin
                currentMin = mat(x, j, z);
                y = j;
                minUpdated = true;
            end
        end

        % Search along the z-dimension (fix x and y)
        for k = 1:dimZ
            if mat(x, y, k) < currentMin
                currentMin = mat(x, y, k);
                z = k;
                minUpdated = true;
            end
        end
    end

    % Return the indices of the minimum value
end
