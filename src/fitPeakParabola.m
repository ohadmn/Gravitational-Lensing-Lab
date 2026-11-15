function r = fitPeakParabola(lc, opts)
%FITPEAKPARABOLA Time and height of the light-curve maximum from a parabola.
%
%   r = fitPeakParabola(lc) fits a weighted parabola to the points around
%   the peak and returns the vertex with propagated 1-sigma errors.
%
%   Options
%     'Window'        [tMin tMax] time window to fit. Default: chosen
%                     automatically as the contiguous run of points around
%                     the peak whose (smoothed) flux exceeds
%                     1 + PeakFraction * (peak - 1).
%     'PeakFraction'  default 0.8
%     'MinPoints'     minimum number of points in the automatic window (7)
%
%   Output struct: tPeak, tPeakErr, fluxPeak, fluxPeakErr, coeffs, cov,
%   chi2, dof, index (points used), tRef (time origin of the coefficients).

    arguments
        lc struct
        opts.Window double = []
        opts.PeakFraction (1,1) double = 0.8
        opts.MinPoints (1,1) double = 7
    end

    if isempty(opts.Window)
        smooth = movmedian(lc.flux, 5);
        [peak, iPeak] = max(smooth);
        frac = opts.PeakFraction;
        while true
            level = 1 + frac * (peak - 1);
            iL = iPeak; while iL > 1 && smooth(iL - 1) > level, iL = iL - 1; end
            iR = iPeak; while iR < numel(smooth) && smooth(iR + 1) > level, iR = iR + 1; end
            if iR - iL + 1 >= opts.MinPoints || frac <= 0.2
                break
            end
            frac = frac - 0.1;   % too few points: widen the window
        end
        index = (iL:iR).';
    else
        index = find(lc.time >= opts.Window(1) & lc.time <= opts.Window(2));
    end
    if numel(index) < 4
        error('fitPeakParabola:tooFewPoints', 'Only %d points in the peak window.', numel(index));
    end

    % Shift the time origin to keep the normal equations well conditioned
    % (raw Julian dates ~2.45e6 raised to the 4th power lose all precision).
    tRef = mean(lc.time(index));
    x = lc.time(index) - tRef;
    [p, C, chi2] = polyfitWeighted(x, lc.flux(index), lc.fluxErr(index), 2);
    a = p(1); b = p(2); c = p(3);
    if a >= 0
        warning('fitPeakParabola:notAMaximum', 'Fitted parabola opens upwards; no maximum in the window.');
    end

    xPeak = -b / (2 * a);
    Jt = [b / (2 * a^2), -1 / (2 * a), 0];        % d(xPeak)/d(a,b,c)
    Jf = [b^2 / (4 * a^2), -b / (2 * a), 1];      % d(fPeak)/d(a,b,c), fPeak = c - b^2/(4a)

    r.tPeak = xPeak + tRef;
    r.tPeakErr = sqrt(Jt * C * Jt.');
    r.fluxPeak = c - b^2 / (4 * a);
    r.fluxPeakErr = sqrt(Jf * C * Jf.');
    r.coeffs = p;
    r.cov = C;
    r.chi2 = chi2;
    r.dof = numel(index) - 3;
    r.index = index;
    r.tRef = tRef;
end
