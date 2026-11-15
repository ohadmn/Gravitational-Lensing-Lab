function u = uFromMagnification(mu)
%UFROMMAGNIFICATION Invert mu(u) = (u^2+2)/(u sqrt(u^2+4)) for u (mu > 1).
    mu = max(mu, 1 + 1e-6);
    u = sqrt(2 * (mu ./ sqrt(mu.^2 - 1) - 1));
end
