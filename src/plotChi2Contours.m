function fig = plotChi2Contours(res, opts)
%PLOTCHI2CONTOURS Delta-chi2 confidence contours for every pair of fitted parameters.
%
%   fig = plotChi2Contours(res) takes the output of gridFit and draws one
%   panel per pair of free parameters (1 panel for 2 parameters, 3 for 3,
%   6 for 4).
%
%   Options
%     'Mode'    'profile' (default): chi2 is minimised over the parameters
%               not shown, which gives the true joint confidence regions.
%               'slice': the other parameters are held at their best-fit
%               values (what the original scripts drew; regions are
%               smaller because correlations are ignored).
%     'Levels'  delta-chi2 levels (default [1 2.3 2.71 6.63]:
%               1    -> projections give the 1-sigma error of one parameter
%               2.3  -> 68.3% joint region for two parameters
%               2.71 -> 90% for one parameter
%               6.63 -> 99% for one parameter)
%     'Title'   super-title

    arguments
        res struct
        opts.Mode (1,1) string {mustBeMember(opts.Mode, ["profile" "slice"])} = "profile"
        opts.Levels double = [1 2.3 2.71 6.63]
        opts.Title (1,1) string = ""
    end

    names = psplParamNames();
    labels = struct('T0', 'T_0', 'tau', '\tau [days]', 'uMin', 'u_{min}', 'fbl', 'f_{bl}');
    dims = arrayfun(@(f) find(names == f), res.free);
    if numel(dims) < 2
        error('plotChi2Contours:needTwoParams', 'Contours need at least two free parameters.');
    end
    pairs = nchoosek(dims, 2);

    [~, bestSubs] = gridMinimum(res.chi2);
    bestSubs(end+1:4) = 1;
    delta = res.chi2 - res.chi2min;

    fig = figure('Color', 'w');
    nPairs = size(pairs, 1);
    nCols = min(nPairs, 3);
    tiles = tiledlayout(fig, ceil(nPairs / nCols), nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig.Position(3:4) = [420 * nCols, 360 * ceil(nPairs / nCols)];

    for k = 1:nPairs
        dx = pairs(k, 1);  dy = pairs(k, 2);
        Z = delta;
        for other = setdiff(1:4, [dx dy])
            if opts.Mode == "profile"
                Z = min(Z, [], other);
            else
                index = repmat({':'}, 1, 4);
                index{other} = bestSubs(other);
                Z = Z(index{:});
            end
        end
        Z = squeeze(Z);                       % rows: dx, columns: dy

        [x, xLabel] = axisValues(res, names(dx), labels);
        [y, yLabel] = axisValues(res, names(dy), labels);

        ax = nexttile(tiles);
        contour(ax, x, y, Z.', opts.Levels, 'ShowText', true, 'LineWidth', 1);
        hold(ax, 'on');
        plot(ax, x(bestSubs(dx)), y(bestSubs(dy)), 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
        hold(ax, 'off');
        xlabel(ax, xLabel);  ylabel(ax, yLabel);
        grid(ax, 'on');
    end
    if opts.Title ~= ""
        title(tiles, opts.Title);
    end
end

function [values, label] = axisValues(res, name, labels)
    values = res.grid.(name);
    label = labels.(name);
    if name == "T0"
        % Show T0 relative to a round reference so the tick labels are readable.
        ref = floor(res.best.T0);
        values = values - ref;
        label = sprintf('T_0 - %d [days]', ref);
    end
end
