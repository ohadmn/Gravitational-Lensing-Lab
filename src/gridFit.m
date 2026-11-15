function res = gridFit(lc, guess, free, opts)
%GRIDFIT Brute-force chi-squared fit of any subset of the PSPL parameters.
%
%   res = gridFit(lc, guess, free) scans a regular grid around the starting
%   values in the struct GUESS (fields T0, tau, uMin and optionally fbl)
%   for the parameters named in FREE, e.g. ["tau" "uMin"] or
%   ["T0" "tau" "uMin" "fbl"]. The remaining parameters stay fixed at their
%   GUESS values.
%
%   The grid is refined automatically: after each pass it is re-centred on
%   the best point and rescaled so that it spans a few sigma of every free
%   parameter (and it is widened when the minimum or the 1-sigma region
%   touches an edge). The result therefore does not depend on hand-tuned
%   search ranges for each event.
%
%   Options
%     'N'          grid points per free parameter (default 300 / 100 / 60
%                  for 1-2 / 3 / 4 free parameters)
%     'RelSpan'    struct with the initial half-width of each axis:
%                  tau, uMin, fbl relative to their guess; T0 in units of tau
%                  (default tau 0.5, uMin 0.5, fbl 0.25, T0 0.15)
%     'MaxRefine'  maximum number of refinement passes (default 5)
%     'SigmaSpan'  half-width of the final grid in units of sigma (default 5)
%     'Verbose'    print progress (default true)
%
%   Output struct
%     best      struct of best-fit parameters (fixed ones included)
%     errLow, errHigh   1-sigma errors of the free parameters, from the
%               range where the profile chi2 (minimised over all other
%               free parameters) rises by 1 above the minimum
%     chi2min, dof
%     free      names of the free parameters
%     grid      struct of axis vectors of the final grid
%     chi2      4-D chi2 array, dimensions ordered [T0, tau, uMin, fbl]
%     edgeHit   true if the 1-sigma region still touches the grid edge
%
%   See also chi2Grid, plotChi2Contours, fitPSPL.

    arguments
        lc struct
        guess struct
        free string
        opts.N double = []
        opts.RelSpan struct = struct()
        opts.MaxRefine (1,1) double = 5
        opts.SigmaSpan (1,1) double = 5
        opts.Verbose (1,1) logical = true
    end

    names = psplParamNames();
    free = intersect(names, free, 'stable');     % canonical order
    if isempty(free)
        error('gridFit:noFreeParams', 'FREE must name at least one of: %s', strjoin(names, ', '));
    end
    if ~isfield(guess, 'fbl'), guess.fbl = 1; end

    if isempty(opts.N)
        defaults = [300 300 100 60];
        opts.N = defaults(numel(free));
    end

    span = struct('T0', 0.15, 'tau', 0.5, 'uMin', 0.5, 'fbl', 0.25);
    for f = string(fieldnames(opts.RelSpan)).'
        span.(f) = opts.RelSpan.(f);
    end

    centre = guess;
    halfWidth = struct('T0', span.T0 * guess.tau, 'tau', span.tau * guess.tau, ...
                       'uMin', span.uMin * guess.uMin, 'fbl', span.fbl * guess.fbl);

    for pass = 1:opts.MaxRefine
        grid = buildGrid(centre, halfWidth, free, opts.N);
        chi2 = chi2Grid(lc, grid);
        [chi2min, subs] = gridMinimum(chi2);
        subs(end+1:4) = 1;

        converged = true;
        edgeHit = false;
        for f = free
            dim = find(names == f);
            values = grid.(f);
            best.(f) = values(subs(dim));
            [lo, hi, atEdge] = profileInterval(chi2, dim, values, chi2min);
            errLow.(f) = best.(f) - lo;
            errHigh.(f) = hi - best.(f);

            if atEdge
                newHalfWidth = 2 * halfWidth.(f);         % region cut off: widen
                edgeHit = true;
            else
                sigma = max([errLow.(f), errHigh.(f), 2 * (values(2) - values(1))]);
                newHalfWidth = opts.SigmaSpan * sigma;    % zoom to a few sigma
            end
            if abs(newHalfWidth / halfWidth.(f) - 1) > 0.3
                converged = false;
            end
            halfWidth.(f) = newHalfWidth;
            centre.(f) = best.(f);
        end
        if opts.Verbose
            fprintf('  grid pass %d (%s): chi2 = %.2f\n', pass, strjoin(free, ', '), chi2min);
        end
        if converged, break, end
    end
    if edgeHit
        warning('gridFit:edgeHit', ...
            'The 1-sigma region still touches the grid edge; errors are lower limits.');
    end

    for f = setdiff(names, free)
        best.(f) = guess.(f);
    end
    res.best = orderfields(best, cellstr(names));
    res.errLow = errLow;
    res.errHigh = errHigh;
    res.chi2min = chi2min;
    res.dof = numel(lc.time) - numel(free);
    res.free = free;
    res.grid = grid;
    res.chi2 = chi2;
    res.edgeHit = edgeHit;
end

function grid = buildGrid(centre, halfWidth, free, n)
    for f = psplParamNames()
        if ismember(f, free)
            lo = centre.(f) - halfWidth.(f);
            hi = centre.(f) + halfWidth.(f);
            if f ~= "T0"
                % tau, uMin and fbl are positive quantities.
                lo = max(lo, 1e-3 * max(centre.(f), eps));
            end
            grid.(f) = linspace(lo, hi, n);
        else
            grid.(f) = centre.(f);
        end
    end
end

function [lo, hi, atEdge] = profileInterval(chi2, dim, values, chi2min)
    % Profile: minimise chi2 over every dimension except DIM.
    profile = chi2;
    for other = setdiff(1:4, dim)
        profile = min(profile, [], other);
    end
    delta = profile(:).' - chi2min;
    inside = find(delta <= 1);
    iLo = inside(1);
    iHi = inside(end);
    atEdge = iLo == 1 || iHi == numel(values);
    lo = crossing(values, delta, iLo, iLo - 1);
    hi = crossing(values, delta, iHi, iHi + 1);
end

function v = crossing(values, delta, iIn, iOut)
    % Linear interpolation of the point where delta chi2 = 1.
    if iOut < 1 || iOut > numel(values)
        v = values(iIn);
    else
        frac = (1 - delta(iIn)) / (delta(iOut) - delta(iIn));
        v = values(iIn) + frac * (values(iOut) - values(iIn));
    end
end
