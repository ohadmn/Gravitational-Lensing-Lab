function out = analyzeEvent(lc, opts)
%ANALYZEEVENT Full microlensing analysis of one light curve.
%
%   out = analyzeEvent(lc) runs, in order:
%     1. parabola fit of the peak            -> T0 and peak flux
%     2. 2-D chi2 grid over (tau, uMin)      -> T0 fixed from step 1
%     3. 3-D chi2 grid over (T0, tau, uMin)
%     4. 4-D chi2 grid over (T0, tau, uMin, fbl)   [only if 'Blended']
%     5. non-linear fit with fminsearch (same free parameters as the
%        largest grid) as an independent cross-check
%     6. lens mass from tau
%   and draws the light curve and the delta-chi2 contours of every grid.
%
%   LC is a struct from loadLightCurve / makeSyntheticLightCurve.
%
%   Options
%     'Name'       label for printing and file names ("event")
%     'Blended'    also fit the blending fraction fbl (false)
%     'Guess'      struct with any of T0, tau, uMin, fbl to start from;
%                  missing fields are estimated from the data
%     'Reference'  struct of published values to compare with (e.g. OGLE)
%     'GridN'      [N2 N3 N4] grid points per axis for the 2/3/4-D grids
%                  (default [300 100 60])
%     'SaveDir'    folder for PNG figures ("" = do not save)
%     'ShowPlots'  create figures (true)
%
%   Output struct: peak, grid2, grid3, grid4 (if blended), simplex, mass.

    arguments
        lc struct
        opts.Name (1,1) string = "event"
        opts.Blended (1,1) logical = false
        opts.Guess struct = struct()
        opts.Reference struct = struct()
        opts.GridN (1,3) double = [300 100 60]
        opts.SaveDir (1,1) string = ""
        opts.ShowPlots (1,1) logical = true
    end

    fprintf('\n===== %s  (%d points, baseline %.3f mag) =====\n', ...
        opts.Name, numel(lc.time), lc.baselineMag);

    % ---- starting values -------------------------------------------------
    guess = estimateInitialParams(lc);
    for f = string(fieldnames(opts.Guess)).'
        guess.(f) = opts.Guess.(f);
    end

    % ---- 1. parabola around the peak --------------------------------------
    out.peak = fitPeakParabola(lc);
    fprintf('Parabola fit of the peak (%d points):\n', numel(out.peak.index));
    fprintf('  T0        = %.3f +- %.3f\n', out.peak.tPeak, out.peak.tPeakErr);
    fprintf('  peak flux = %.3f +- %.3f  (chi2/dof = %.1f/%d)\n', ...
        out.peak.fluxPeak, out.peak.fluxPeakErr, out.peak.chi2, out.peak.dof);
    if ~isfield(opts.Guess, 'T0') && out.peak.coeffs(1) < 0 ...
            && abs(out.peak.tPeak - guess.T0) < guess.tau
        guess.T0 = out.peak.tPeak;
    end

    % ---- 2-4. chi2 grids ---------------------------------------------------
    out.grid2 = gridFit(lc, guess, ["tau" "uMin"], 'N', opts.GridN(1));
    report(sprintf('2-D grid (tau, uMin; T0 = %.3f, fbl = %g fixed)', guess.T0, guess.fbl), out.grid2, opts.Reference);

    out.grid3 = gridFit(lc, out.grid2.best, ["T0" "tau" "uMin"], 'N', opts.GridN(2));
    report(sprintf('3-D grid (T0, tau, uMin; fbl = %g fixed)', guess.fbl), out.grid3, opts.Reference);
    final = out.grid3;

    if opts.Blended
        out.grid4 = gridFit(lc, out.grid3.best, ["T0" "tau" "uMin" "fbl"], 'N', opts.GridN(3));
        report('4-D grid (T0, tau, uMin, fbl)', out.grid4, opts.Reference);
        final = out.grid4;
    end

    % ---- 5. simplex cross-check -------------------------------------------
    out.simplex = fitPSPL(lc, guess, final.free);
    fprintf('fminsearch (%s):  chi2/dof = %.1f/%d\n', ...
        strjoin(final.free, ', '), out.simplex.chi2min, out.simplex.dof);
    for f = out.simplex.free
        fprintf('  %-5s = %s +- %s\n', f, fmt(f, out.simplex.best.(f)), fmt("", out.simplex.err.(f)));
    end

    % ---- 6. lens mass -------------------------------------------------------
    tau = final.best.tau;
    out.mass = lensMass(tau);
    massRange = lensMass([tau - final.errLow.tau, tau + final.errHigh.tau]);
    fprintf('Lens mass (v = 150 km/s, Ds = 8 kpc, x = 1/2): %.3f (+%.3f / -%.3f) solar masses\n', ...
        out.mass, massRange(2) - out.mass, out.mass - massRange(1));

    % ---- figures -------------------------------------------------------------
    if opts.ShowPlots
        figs.lightcurve = plotLightCurve(lc, final.best, 'Title', opts.Name + ": light curve and best fit");
        figs.contours2D = plotChi2Contours(out.grid2, 'Title', opts.Name + ": \Delta\chi^2, (\tau, u_{min}) grid");
        figs.contours3D = plotChi2Contours(out.grid3, 'Title', opts.Name + ": \Delta\chi^2, 3-parameter grid");
        if opts.Blended
            figs.contours4D = plotChi2Contours(out.grid4, 'Title', opts.Name + ": \Delta\chi^2, 4-parameter grid");
        end
        if opts.SaveDir ~= ""
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            for f = string(fieldnames(figs)).'
                exportgraphics(figs.(f), fullfile(opts.SaveDir, opts.Name + "_" + f + ".png"), 'Resolution', 150);
            end
        end
    end
end

function report(label, res, reference)
    fprintf('%s:  chi2/dof = %.1f/%d\n', label, res.chi2min, res.dof);
    for f = res.free
        line = sprintf('  %-5s = %s  (+%s / -%s)', f, fmt(f, res.best.(f)), ...
            fmt("", res.errHigh.(f)), fmt("", res.errLow.(f)));
        if isfield(reference, f)
            sigma = max(res.errHigh.(f), res.errLow.(f));
            line = sprintf('%s   reference %s  (%.1f sigma away)', line, ...
                fmt(f, reference.(f)), abs(res.best.(f) - reference.(f)) / sigma);
        end
        fprintf('%s\n', line);
    end
end

function s = fmt(name, value)
    if name == "T0"
        s = sprintf('%.3f', value);
    else
        s = sprintf('%.4g', value);
    end
end
