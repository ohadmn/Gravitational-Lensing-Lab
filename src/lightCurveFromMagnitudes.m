function lc = lightCurveFromMagnitudes(time, mag, magErr, baselineMag)
%LIGHTCURVEFROMMAGNITUDES Build a normalised-flux light-curve struct.
%
%   lc = lightCurveFromMagnitudes(time, mag, magErr) estimates the baseline
%   magnitude from the data; pass baselineMag as a 4th argument to fix it.
%
%   Conversion (F0 cancels out in the normalisation):
%       F / F_base = 10^(-0.4 (m - m_base))
%       sigma_F    = 0.4 ln(10) * sigma_m * F
%
%   See also loadLightCurve.

    if nargin < 4, baselineMag = NaN; end

    time = time(:); mag = mag(:); magErr = magErr(:);
    ok = isfinite(time) & isfinite(mag) & isfinite(magErr) & magErr > 0;
    [time, order] = sort(time(ok));
    mag = mag(ok);       mag = mag(order);
    magErr = magErr(ok); magErr = magErr(order);
    if numel(time) < 5
        error('lightCurveFromMagnitudes:tooFewPoints', 'Fewer than 5 valid data points.');
    end

    lc = struct('time', time, 'mag', mag, 'magErr', magErr);

    if isnan(baselineMag)
        % Pass 1: crude baseline = median of everything.
        lc = normalise(lc, median(mag));
        % Pass 2: locate the event and use only points far away from it.
        guess = estimateInitialParams(lc);
        far = abs(time - guess.T0) > 5 * guess.tau;
        if nnz(far) >= 10
            baselineMag = median(mag(far));
        else
            warning('lightCurveFromMagnitudes:noBaseline', ...
                ['Fewer than 10 points lie outside the event; the baseline is ', ...
                 'the median of all points and may be biased. Consider passing BaselineMag.']);
            baselineMag = median(mag);
        end
    end
    lc = normalise(lc, baselineMag);
end

function lc = normalise(lc, baselineMag)
    lc.baselineMag = baselineMag;
    lc.flux = 10.^(-0.4 * (lc.mag - baselineMag));
    lc.fluxErr = 0.4 * log(10) * lc.magErr .* lc.flux;
end
