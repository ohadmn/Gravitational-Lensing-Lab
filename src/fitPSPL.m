function res = fitPSPL(lc, guess, free)
%FITPSPL Non-linear chi-squared fit of the PSPL model with fminsearch.
%
%   res = fitPSPL(lc, guess) fits T0, tau, uMin and fbl starting from the
%   struct GUESS. res = fitPSPL(lc, guess, free) fits only the parameters
%   named in FREE (e.g. ["T0" "tau" "uMin"] for an unblended event).
%
%   Internally T0 is fitted as an offset from the starting value, measured
%   in units of tau, so that all simplex coordinates are of order one
%   (a simplex on raw Julian dates ~2.45e6 does not converge).
%
%   Output struct: best, err (1-sigma, from the curvature of chi2 at the
%   minimum: cov = 2 inv(Hessian)), chi2min, dof, free, cov.
%
%   See also gridFit.

    arguments
        lc struct
        guess struct
        free string = psplParamNames()
    end

    names = psplParamNames();
    free = intersect(names, free, 'stable');
    if ~isfield(guess, 'fbl'), guess.fbl = 1; end

    % Scaled coordinates: x = (p - offset) / scale.
    offset = struct('T0', guess.T0, 'tau', 0, 'uMin', 0, 'fbl', 0);
    scale  = struct('T0', guess.tau, 'tau', guess.tau, 'uMin', guess.uMin, 'fbl', guess.fbl);

    x0 = zeros(1, numel(free));
    for i = 1:numel(free)
        x0(i) = (guess.(free(i)) - offset.(free(i))) / scale.(free(i));
    end

    options = optimset('TolX', 1e-9, 'TolFun', 1e-9, 'MaxIter', 1e5, 'MaxFunEvals', 1e5);
    x = x0;
    for restart = 1:3            % restarts guard against a stalled simplex
        [x, chi2min] = fminsearch(@objective, x, options);
    end

    best = unpack(x);
    best.tau = abs(best.tau);    % the model only depends on |tau|, |uMin|
    best.uMin = abs(best.uMin);

    % Covariance from a central-difference Hessian in the scaled coordinates.
    n = numel(x);
    h = 1e-3 * max(abs(x), 0.1);
    H = zeros(n);
    for i = 1:n
        for j = i:n
            ei = zeros(1, n); ei(i) = h(i);
            ej = zeros(1, n); ej(j) = h(j);
            H(i, j) = (objective(x + ei + ej) - objective(x + ei - ej) ...
                     - objective(x - ei + ej) + objective(x - ei - ej)) / (4 * h(i) * h(j));
            H(j, i) = H(i, j);
        end
    end
    s = arrayfun(@(f) scale.(f), free);
    cov = 2 * inv(H) .* (s.' * s); %#ok<MINV>
    for i = 1:n
        if all(isfinite(cov(:))) && cov(i, i) > 0
            err.(free(i)) = sqrt(cov(i, i));
        else
            err.(free(i)) = NaN;
        end
    end

    res.best = best;
    res.err = err;
    res.chi2min = chi2min;
    res.dof = numel(lc.time) - n;
    res.free = free;
    res.cov = cov;

    function p = unpack(x)
        p = guess;
        for k = 1:numel(free)
            p.(free(k)) = offset.(free(k)) + scale.(free(k)) * x(k);
        end
    end

    function c = objective(x)
        c = sum(((lc.flux - psplMagnification(lc.time, unpack(x))) ./ lc.fluxErr).^2);
    end
end
