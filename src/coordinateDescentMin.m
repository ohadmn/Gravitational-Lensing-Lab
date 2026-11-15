function [minValue, subs] = coordinateDescentMin(A, start)
%COORDINATEDESCENTMIN Minimum of an N-D array by axis-by-axis line searches.
%
%   [minValue, subs] = coordinateDescentMin(A) starts at the first element
%   and repeatedly moves to the smallest element along each axis in turn,
%   until a full sweep over all axes brings no improvement.
%   coordinateDescentMin(A, start) starts from the subscripts START.
%
%   This is the N-dimensional generalisation of the lab's original
%   matrixMinimization_2D/3D/4D routines. It only inspects a few lines of
%   the array, but it converges to a point that is minimal along every
%   axis - which is the global minimum only for well-behaved (e.g. convex,
%   single-valley) surfaces. Use gridMinimum when the whole array is
%   available anyway.

    nd = ndims(A);
    if nargin < 2, start = ones(1, nd); end
    subs = num2cell(start);
    minValue = A(subs{:});

    improved = true;
    while improved
        improved = false;
        for dim = 1:nd
            line = subs;
            line{dim} = ':';
            [value, index] = min(A(line{:}), [], 'all');
            if value < minValue
                minValue = value;
                subs{dim} = index;
                improved = true;
            end
        end
    end
    subs = cell2mat(subs);
end
