function [p, C, chi2] = polyfitWeighted(x, y, sigma, degree)
%POLYFITWEIGHTED Weighted linear least-squares polynomial fit.
%
%   [p, C, chi2] = polyfitWeighted(x, y, sigma, degree) minimises
%       chi2 = sum(((y - polyval(p, x)) ./ sigma).^2)
%   and returns the coefficients p (highest power first, as in POLYVAL),
%   their covariance matrix C and the minimum chi2.
%
%   degree = 2 reproduces the parabola fit of the lab's first part.

    x = x(:); y = y(:); sigma = sigma(:);
    if numel(x) <= degree
        error('polyfitWeighted:tooFewPoints', ...
            'Need more than %d points for a degree-%d fit.', degree, degree);
    end
    A = (x .^ (degree:-1:0)) ./ sigma;   % weighted design matrix
    b = y ./ sigma;
    p = (A \ b).';
    C = inv(A.' * A);
    chi2 = sum((b - A * p.').^2);
end
