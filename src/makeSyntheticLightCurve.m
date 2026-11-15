function lc = makeSyntheticLightCurve(p, opts)
%MAKESYNTHETICLIGHTCURVE Simulated OGLE-like light curve with known truth.
%
%   lc = makeSyntheticLightCurve(p) draws noisy magnitudes from the PSPL
%   model with parameter struct p (T0, tau, uMin, optional fbl). Useful for
%   testing the fitting code and for running the demo without data files.
%
%   Options: 'NumPoints' (600), 'Span' total duration in days (default
%   covers T0 +- 15 tau, at least 400 d), 'BaselineMag' (18), 'MagErr'
%   typical photometric error at baseline (0.02), 'Seed' (1).

    arguments
        p struct
        opts.NumPoints (1,1) double = 600
        opts.Span (1,1) double = NaN
        opts.BaselineMag (1,1) double = 18
        opts.MagErr (1,1) double = 0.02
        opts.Seed (1,1) double = 1
    end

    stream = RandStream('mt19937ar', 'Seed', opts.Seed);
    if isnan(opts.Span), opts.Span = max(400, 30 * p.tau); end

    % Irregular sampling, like ground-based survey data.
    time = sort(p.T0 + opts.Span * (rand(stream, opts.NumPoints, 1) - 0.5));
    flux = psplMagnification(time, p);
    magErr = opts.MagErr ./ sqrt(flux);                  % brighter = more precise
    mag = opts.BaselineMag - 2.5 * log10(flux) + magErr .* randn(stream, opts.NumPoints, 1);

    lc = lightCurveFromMagnitudes(time, mag, magErr);
    lc.truth = p;
end
