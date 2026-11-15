function lc = loadLightCurve(filePath, opts)
%LOADLIGHTCURVE Read a photometry file and convert magnitudes to normalised flux.
%
%   lc = loadLightCurve(filePath) reads a whitespace-delimited text file whose
%   first three columns are time [days, e.g. HJD], magnitude and magnitude
%   error (the OGLE EWS "phot.dat" layout; extra columns are ignored).
%
%   The flux is normalised to the un-lensed baseline, so that far from the
%   event flux = 1 and during the event flux = 1 - fbl + fbl * mu(t).
%
%   lc = loadLightCurve(..., 'BaselineMag', m0) uses a known baseline
%   magnitude. Otherwise the baseline is the median magnitude of all points
%   further than 5 Einstein-crossing times from the peak (found from a
%   first-pass estimate of the event), which is far more robust than
%   normalising by a single data point.
%
%   Output struct fields:
%     time, mag, magErr   raw columns (sorted by time, invalid rows removed)
%     flux, fluxErr       baseline-normalised flux and its 1-sigma error
%     baselineMag         baseline magnitude used for the normalisation
%
%   See also makeSyntheticLightCurve, estimateInitialParams.

    arguments
        filePath (1,1) string
        opts.BaselineMag (1,1) double = NaN
    end

    if ~isfile(filePath)
        error('loadLightCurve:fileNotFound', 'Data file not found: %s', filePath);
    end
    data = readmatrix(filePath, 'FileType', 'text');
    if size(data, 2) < 3
        error('loadLightCurve:badFormat', ...
            'Expected at least 3 columns (time, magnitude, error) in %s', filePath);
    end
    lc = lightCurveFromMagnitudes(data(:,1), data(:,2), data(:,3), opts.BaselineMag);
end
