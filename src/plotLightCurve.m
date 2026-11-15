function fig = plotLightCurve(lc, p, opts)
%PLOTLIGHTCURVE Data with error bars, and optionally a model curve.
%
%   fig = plotLightCurve(lc) plots the normalised flux.
%   fig = plotLightCurve(lc, p) overlays the PSPL model with parameters p
%   and zooms on the event (T0 +- 4 tau). Option 'Title'.

    arguments
        lc struct
        p struct = struct([])
        opts.Title (1,1) string = ""
    end

    % Plot relative to a round reference date so the tick labels are readable.
    if ~isempty(p), tRef = floor(p.T0); else, tRef = floor(lc.time(1)); end

    fig = figure('Color', 'w');
    errorbar(lc.time - tRef, lc.flux, lc.fluxErr, '.', 'MarkerSize', 8, 'CapSize', 0, ...
        'Color', [0.35 0.35 0.35], 'DisplayName', 'data');
    hold on
    if ~isempty(p)
        tModel = linspace(p.T0 - 4 * p.tau, p.T0 + 4 * p.tau, 1000);
        plot(tModel - tRef, psplMagnification(tModel, p), 'r-', 'LineWidth', 1.5, 'DisplayName', 'PSPL model');
        xlim([tModel(1), tModel(end)] - tRef);
        legend('Location', 'northeast');
    end
    hold off
    xlabel(sprintf('Time - %d [days]', tRef));
    ylabel('Flux / baseline flux');
    grid on
    if opts.Title ~= "", title(opts.Title); end
end
