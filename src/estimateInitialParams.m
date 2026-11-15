function p = estimateInitialParams(lc)
%ESTIMATEINITIALPARAMS Rough (T0, tau, uMin) directly from a light curve.
%
%   p = estimateInitialParams(lc) needs no prior knowledge of the event, so
%   the analysis can run on any light curve without hard-coded catalogue
%   values. Assumes no blending (p.fbl = 1).
%
%     T0    time of the maximum of the median-smoothed flux
%     uMin  from the peak magnification
%     tau   from the half-maximum width of the peak
%
%   The values are only starting points for gridFit / fitPSPL.

    smooth = movmedian(lc.flux, 5);
    [peak, iPeak] = max(smooth);
    peak = max(peak, 1.02);

    p.T0 = lc.time(iPeak);
    p.uMin = uFromMagnification(peak);

    % Walk outwards from the peak while the curve is above half maximum.
    half = 1 + (peak - 1) / 2;
    iL = iPeak; while iL > 1 && smooth(iL - 1) > half, iL = iL - 1; end
    iR = iPeak; while iR < numel(smooth) && smooth(iR + 1) > half, iR = iR + 1; end
    halfWidth = (lc.time(iR) - lc.time(iL)) / 2;
    if halfWidth <= 0
        halfWidth = median(diff(lc.time));
    end
    uHalf = uFromMagnification(half);
    p.tau = halfWidth / sqrt(uHalf^2 - p.uMin^2);
    p.fbl = 1;
end
