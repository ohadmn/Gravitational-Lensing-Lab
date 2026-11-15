function chi2 = chi2Grid(lc, grid)
%CHI2GRID Chi-squared of the PSPL model on an N-dimensional parameter grid.
%
%   chi2 = chi2Grid(lc, grid) evaluates
%       chi2 = sum_k ((flux_k - model(t_k)) / sigma_k)^2
%   for every combination of the values in the struct GRID, whose fields
%   T0, tau, uMin and fbl are each a vector (a scalar keeps that parameter
%   fixed). The result is a 4-D array of size
%       [numel(T0), numel(tau), numel(uMin), numel(fbl)]
%   so any 1- to 4-parameter search uses the same function.
%
%   The model is linear in fbl:  model - 1 = fbl * (mu - 1).  Therefore
%       chi2 = Sdd - 2 fbl Sdm + fbl^2 Smm
%   where the sums S only depend on (T0, tau, uMin). The data loop runs
%   over a 3-D grid at most and the fbl axis is added at no extra cost.
%
%   See also gridFit, gridMinimum.

    T0  = reshape(grid.T0,  [], 1);
    tau = reshape(grid.tau,  1, []);
    u2  = reshape(grid.uMin, 1, 1, []).^2;
    fbl = reshape(grid.fbl,  1, 1, 1, []);

    w = 1 ./ lc.fluxErr.^2;
    d = lc.flux - 1;

    Sdd = sum(w .* d.^2);
    Sdm = zeros(numel(T0), numel(tau), numel(u2));
    Smm = Sdm;
    for k = 1:numel(lc.time)
        uSq = u2 + ((lc.time(k) - T0) ./ tau).^2;
        m = (uSq + 2) ./ sqrt(uSq .* (uSq + 4)) - 1;   % mu - 1
        Sdm = Sdm + (w(k) * d(k)) * m;
        Smm = Smm + w(k) * m.^2;
    end
    chi2 = Sdd - 2 * fbl .* Sdm + fbl.^2 .* Smm;
end
