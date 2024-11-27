function [w, x, y, z] = matrixMinimization_4D(mat)
    % Initialize variables
    [dimW, dimX, dimY, dimZ] = size(mat);
    currentMin = flintmax; % Initialize to the largest possible value
    w = 1; x = 1; y = 1; z = 1;
    minUpdated = true;

    % Loop until no smaller value is found
    while minUpdated
        minUpdated = false;

        % Search along the w-dimension fix:(x, y, z)
        for i = 1:dimW
            if mat(i, x, y, z) < currentMin
                currentMin = mat(i, x, y, z);
                w = i;
                minUpdated = true;
            end
        end

        % Search along the x-dimension fix:(w, y, z)
        for j = 1:dimX
            if mat(w, j, y, z) < currentMin
                currentMin = mat(w, j, y, z);
                x = j;
                minUpdated = true;
            end
        end

        % Search along the y-dimension fix:(w, x, z)
        for k = 1:dimY
            if mat(w, x, k, z) < currentMin
                currentMin = mat(w, x, k, z);
                y = k;
                minUpdated = true;
            end
        end

        % Search along the z-dimension fix:(w, x, y)
        for l = 1:dimZ
            if mat(w, x, y, l) < currentMin
                currentMin = mat(w, x, y, l);
                z = l;
                minUpdated = true;
            end
        end
    end

    % Return the indices of the minimum value
end
