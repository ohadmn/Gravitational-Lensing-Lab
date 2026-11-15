function [massSolar, einsteinRadiusAU] = lensMass(tau, opts)
%LENSMASS Lens mass from the Einstein-radius crossing time.
%
%   massSolar = lensMass(tau) with tau in days. The Einstein radius is
%   R_E = v * tau and
%       R_E^2 = (4 G M / c^2) * Ds * x * (1 - x),   x = Dl / Ds
%   so  M = (v tau)^2 c^2 / (4 G Ds x (1 - x)).
%
%   Options (defaults are the values suggested in the lab manual)
%     'Velocity'        transverse lens velocity [km/s]   (150)
%     'SourceDistance'  Ds [kpc]                          (8)
%     'DistanceRatio'   x = Dl / Ds                       (0.5)
%
%   tau may be an array. The second output is R_E in astronomical units.

    arguments
        tau double
        opts.Velocity (1,1) double = 150
        opts.SourceDistance (1,1) double = 8
        opts.DistanceRatio (1,1) double = 0.5
    end

    G = 6.674e-11;  c = 2.998e8;  solarMass = 1.989e30;   % SI
    kpc = 3.0857e19;  au = 1.496e11;  day = 86400;

    RE = opts.Velocity * 1e3 * tau * day;
    x = opts.DistanceRatio;
    mass = RE.^2 * c^2 / (4 * G * opts.SourceDistance * kpc * x * (1 - x));
    massSolar = mass / solarMass;
    einsteinRadiusAU = RE / au;
end
