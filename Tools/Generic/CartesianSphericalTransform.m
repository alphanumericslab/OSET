function y = CartesianSphericalTransform(x, mode)
% Y = CartesianSphericalTransform(X, mode)
% Transforms between Cartesian and n-dimensional Spherical coordinates
% See: https://en.wikipedia.org/wiki/N-sphere
%
% Inputs:
%   X: m x n matrix of row-wise Cartesian or Spherical coordinates
%       if mode = 'CARTtoSPHR' each row of X is of the form x = [x_1, x_2, ..., x_n]
%       if mode = 'SPHRtoCART' each row of X is of the form x = [r, phi_1, ..., phi_{n-1}]
%   mode:
%       'CARTtoSPHR': Cartesian to Spherical coordinates
%       'SPHRtoCART': Spherical to Cartesian coordinates
% Output:
%   Y: m x n matrix of row-wise Spherical or Cartesian coordinates
%       if mode = 'CARTtoSPHR' each row of Y is of the form y = [r, phi_1, ..., phi_{n-1}]
%       if mode = 'SPHRtoCART' each row of Y is of the form y = [x_1, x_2, ..., x_n]
%
% Open Source Electrophysiological Toolbox, version 3.14
% Copyright (C) 2021  Reza Sameni
% reza.sameni@gmail.com

% Correct if x is a single row-wise vector
if(iscolumn(x))
    x = x';
end

method = 'ACOS_FORMULATION'; % 'ACOT_FORMULATION'. TODO: The ACOT algorithm has sign issues.
if(isequal(mode, 'SPHRtoCART'))
    y = zeros(size(x));
    r = x(:, 1);
    phi = x(:, 2 : end);
    factor = r;
    for k = 1 : size(phi, 2)
        y(:, k) = factor .* cos(phi(:, k));
        factor = factor .* sin(phi(:, k));
    end
    y(:, end) = factor;
elseif(isequal(mode, 'CARTtoSPHR'))
    y = zeros(size(x));
    y(:, 1) = sqrt(sum(x .^ 2, 2));
    if(isequal(method, 'ACOT_FORMULATION'))
        for k = 1 : size(x, 2) - 2
            y(:, k + 1) = acot(x(:, k) ./ sqrt(sum(x(:, k + 1 : end).^2)));
        end
        y(:, end) = 2 * acot(( x(:, end - 1) + sqrt(sum(x(:, end - 1 : end).^2, 2)) ) ./ x(:, end));
    elseif(isequal(method, 'ACOS_FORMULATION'))
        for k = 1 : size(x, 2) - 1
            y(:, k + 1) = acos(x(:, k) ./ sqrt(sum(x(:, k : end).^2, 2)));
        end
        sgn = x(:, end) < 0;
        y(sgn, end) = - y(sgn, end);

    end
elseif(isequal(mode, 'CARTtoSPHR-ISO'))
    if size(x, 2) == 3 % CARTtoSPHR-ISO mode is only available for 3-vectors
        y = zeros(size(x));
        y(:, 1) = sqrt(sum(x .^ 2, 2)); %
        y(:, 2) = acos(x(:, 3) ./ y(:, 1)); % theta
        for m = 1 : size(x, 1) % phi
            if x(m, 1) ~= 0
                y(m, 3) = atan2(x(m, 2), x(m, 1));
            else
                if x(m, 2) > 0
                    y(m, 3) = pi/2;
                elseif x(m, 2) < 0
                    y(m, 3) = -pi/2;
                else
                    y(m, 3) = nan;
                end
            end
        end
    else
        error('CARTtoSPHR-ISO mode is only available for 3-vectors')
    end
else
    error('Undefined transformation mode.');
end

% Correct y dim if x was a single row-wise vector
if(iscolumn(x))
    y = y';
end

