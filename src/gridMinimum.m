function [minValue, subs] = gridMinimum(A)
%GRIDMINIMUM Global minimum of an N-dimensional array and its subscripts.
%
%   [minValue, subs] = gridMinimum(A) returns the smallest element of A and
%   a row vector SUBS with one subscript per dimension of A, so that
%   A(subs(1), subs(2), ...) == minValue. Works for any number of
%   dimensions (replaces the separate 2-D / 3-D / 4-D routines).
%
%   See also coordinateDescentMin.

    [minValue, linearIndex] = min(A(:));
    subs = cell(1, ndims(A));
    [subs{:}] = ind2sub(size(A), linearIndex);
    subs = cell2mat(subs);
end
