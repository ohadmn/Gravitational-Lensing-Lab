function f = psplMagnification(t, p)
%PSPLMAGNIFICATION Point-source point-lens (Paczynski) light curve.
%
%   f = psplMagnification(t, p) returns the baseline-normalised flux at the
%   times t for the parameter struct p with fields
%       T0    time of closest approach
%       tau   Einstein-radius crossing time (same units as t)
%       uMin  impact parameter in Einstein radii
%       fbl   (optional, default 1) fraction of the baseline flux that
%             belongs to the lensed source ("blending" parameter)
%
%       u(t)  = sqrt(uMin^2 + ((t - T0)/tau)^2)
%       mu(u) = (u^2 + 2) / (u sqrt(u^2 + 4))
%       f(t)  = 1 - fbl + fbl * mu
%
%   All fields broadcast against t, so arrays of parameters are allowed.

    if isfield(p, 'fbl'), fbl = p.fbl; else, fbl = 1; end
    u2 = p.uMin.^2 + ((t - p.T0) ./ p.tau).^2;
    mu = (u2 + 2) ./ sqrt(u2 .* (u2 + 4));
    f = 1 - fbl + fbl .* mu;
end
