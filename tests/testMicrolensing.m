function tests = testMicrolensing
%TESTMICROLENSING Unit tests.  Run with:  runtests('tests')
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
end

function testPolyfitRecoversCoefficients(t)
    rng(1);
    truth = [1 2 3];
    x = linspace(1, 5, 100);
    sigma = 0.1 * ones(size(x));
    nTrials = 2000;
    p = zeros(nTrials, 3);  chi2 = zeros(nTrials, 1);  C = [];
    for i = 1:nTrials
        y = polyval(truth, x) + sigma .* randn(size(x));
        [p(i, :), C, chi2(i)] = polyfitWeighted(x, y, sigma, 2);
    end
    % Unbiased, scatter matches the covariance matrix, chi2 ~ dof = 97.
    verifyEqual(t, mean(p), truth, 'AbsTol', 5 * sqrt(diag(C)).' / sqrt(nTrials));
    verifyEqual(t, std(p), sqrt(diag(C)).', 'RelTol', 0.1);
    verifyEqual(t, mean(chi2), 97, 'RelTol', 0.05);
end

function testPolyfitMatchesPolyfit(t)
    x = linspace(-3, 3, 30);
    y = 0.5 * x.^3 - x + 2;
    p = polyfitWeighted(x, y, ones(size(x)), 3);
    verifyEqual(t, p, [0.5 0 -1 2], 'AbsTol', 1e-10);
end

function testGridMinimumAnyDimension(t)
    rng(2);
    for nd = 2:4
        A = rand([7 6 5 4] .* (1:4 <= nd) + (1:4 > nd));
        [value, subs] = gridMinimum(A);
        index = num2cell(subs);
        verifyEqual(t, value, min(A(:)));
        verifyEqual(t, A(index{:}), value);
    end
end

function testCoordinateDescentOnConvexSurface(t)
    [a, b, c, d] = ndgrid(1:20, 1:25, 1:15, 1:10);
    A = (a - 13).^2 + 2 * (b - 4).^2 + (c - 9).^2 + 3 * (d - 7).^2 + 0.5 * (a - 13) .* (b - 4);
    [value, subs] = coordinateDescentMin(A);
    verifyEqual(t, subs, [13 4 9 7]);
    verifyEqual(t, value, 0);
end

function testMagnificationLimits(t)
    p = struct('T0', 0, 'tau', 10, 'uMin', 0.1);
    verifyEqual(t, psplMagnification(0, p), (0.01 + 2) / (0.1 * sqrt(4.01)), 'RelTol', 1e-12);
    verifyEqual(t, psplMagnification(1e4, p), 1, 'AbsTol', 1e-6);
    verifyEqual(t, uFromMagnification(psplMagnification(0, p)), 0.1, 'RelTol', 1e-9);
    p.fbl = 0.3;      % blending dilutes the magnification
    verifyEqual(t, psplMagnification(0, p), 1 + 0.3 * ((0.01 + 2) / (0.1 * sqrt(4.01)) - 1), 'RelTol', 1e-12);
end

function testChi2GridMatchesDirectSum(t)
    truth = struct('T0', 100, 'tau', 20, 'uMin', 0.3, 'fbl', 0.7);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 200);
    grid = struct('T0', [99 100.5], 'tau', [18 20 23], 'uMin', [0.25 0.3], 'fbl', [0.6 0.7 0.8 1]);
    chi2 = chi2Grid(lc, grid);
    verifySize(t, chi2, [2 3 2 4]);
    p = struct('T0', 100.5, 'tau', 23, 'uMin', 0.25, 'fbl', 0.8);
    direct = sum(((lc.flux - psplMagnification(lc.time, p)) ./ lc.fluxErr).^2);
    verifyEqual(t, chi2(2, 3, 1, 3), direct, 'RelTol', 1e-9);
end

function testParabolaFindsPeak(t)
    truth = struct('T0', 2454500.25, 'tau', 30, 'uMin', 0.3);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 1500, 'MagErr', 0.01);
    r = fitPeakParabola(lc);
    verifyLessThan(t, abs(r.tPeak - truth.T0), max(5 * r.tPeakErr, 0.5));
end

function testGridFitRecoversUnblendedEvent(t)
    truth = struct('T0', 2454476.1, 'tau', 51, 'uMin', 0.12, 'fbl', 1);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 600, 'Span', 1500);
    guess = estimateInitialParams(lc);
    res = gridFit(lc, guess, ["T0" "tau" "uMin"], 'N', 60, 'Verbose', false);
    verifyWithinSigma(t, res, truth, 4);
    verifyEqual(t, res.chi2min / res.dof, 1, 'AbsTol', 0.25);
    verifyFalse(t, res.edgeHit);
end

function testGridFitRecoversBlendedEvent(t)
    truth = struct('T0', 2454517.5, 'tau', 17, 'uMin', 0.4, 'fbl', 0.5);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 800, 'MagErr', 0.01);
    guess = estimateInitialParams(lc);
    res = gridFit(lc, guess, ["T0" "tau" "uMin" "fbl"], 'N', 40, 'Verbose', false);
    verifyWithinSigma(t, res, truth, 4);
end

function testSimplexAgreesWithGrid(t)
    truth = struct('T0', 2454517.5, 'tau', 17, 'uMin', 0.4, 'fbl', 0.5);
    lc = makeSyntheticLightCurve(truth, 'NumPoints', 800, 'MagErr', 0.01);
    guess = estimateInitialParams(lc);
    grid = gridFit(lc, guess, ["T0" "tau" "uMin" "fbl"], 'N', 40, 'Verbose', false);
    simplex = fitPSPL(lc, guess);
    verifyLessThanOrEqual(t, simplex.chi2min, grid.chi2min + 1e-6);   % continuous optimum
    verifyLessThan(t, grid.chi2min - simplex.chi2min, 1);             % grid is close to it
    for f = grid.free
        verifyLessThan(t, abs(simplex.best.(f) - grid.best.(f)), ...
            max(grid.errLow.(f), grid.errHigh.(f)));
    end
end

function testLensMass(t)
    % tau = 20 d at 150 km/s: R_E = 1.73 AU, M = 0.18 solar masses.
    [m, re] = lensMass(20);
    verifyEqual(t, re, 1.733, 'RelTol', 1e-3);
    verifyEqual(t, m, 0.184, 'RelTol', 0.01);
    verifyEqual(t, lensMass(40) / m, 4, 'RelTol', 1e-12);     % M ~ tau^2
end

function verifyWithinSigma(t, res, truth, nSigma)
    for f = res.free
        sigma = max(res.errLow.(f), res.errHigh.(f));
        verifyLessThan(t, abs(res.best.(f) - truth.(f)), nSigma * sigma, ...
            sprintf('%s is more than %d sigma from the truth', f, nSigma));
    end
end
