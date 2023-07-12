function [y, B, conv] = easi_source_separation(x, nsou, lambda, nlintype)
% easi_source_separation - An adaptive blind source separation algorithm using the EASI (equivariant adaptive separation by independence) algorithm
%
% Syntax: [y, B, conv] = easi_source_separation(x, nsou, lambda, nlintype)
%
% Inputs:
%   x: input signal as a matrix of size (ncap x T), where ncap is the number of input channels and T is the length of the signal
%   nsou: number of sources
%   lambda: adaptation parameter
%   nlintype: nonlinearity type for the adaptation (y2 denotes y.^2)
%             1: g = diag(y2) .* y(:, t)
%             2: g = y(:, t) .* diag(0.1 * y2) .* diag(y2)
%             3: g = y(:, t) .* sqrt(diag(y2))
%             4: g = y(:, t) .* log(diag(y2))
%             5: g = -y(:, t) .* (diag(y2) < 0.9)
%             6: g = y(:, t) ./ log(diag(y2))
%             7: g = -y(:, t) ./ sqrt(diag(y2))
%             8: g = -y(:, t) ./ diag(y2)
%
% Outputs:
%   y: separated sources as a matrix of size (nsou x T)
%   B: separation matrix of size (nsou x ncap)
%   conv: convergence of the algorithm at each iteration (1 x T)
%
% References:
%   Beate Laheld and Jean-François Cardoso, "Adaptive source separation without prewhitening," Proc. EUSIPCO'94, 183-186, Edinburgh, Sep. 1994.
%   This code has been adapted from a demo file obtained from: http://bsp.teithe.gr/members/downloads/easidemo/easi_demo.m
% 
%   Revision History:
%       2019: First release
%       2023: Added help and renamed from deprecated version EASI()
% 
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

ncap = size(x, 1); % Number of input channels
T = size(x, 2); % Length of the signal
idsou = eye(nsou); % Identity matrix for nsou sources

B = randn(nsou, ncap); % Initialization of the separation matrix
y = zeros(nsou, T); % Initialization of the separated sources
conv = zeros(1, T); % Initialization of the convergence vector

for t = 1:T
    y(:, t) = B * x(:, t); % Separation of sources using the separation matrix
    y2 = y(:, t) * y(:, t)'; % Square of the separated sources
    
    % Nonlinearity adaptation based on the specified nlintype
    switch nlintype
        case 1
            g = diag(y2) .* y(:, t);
        case 2
            g = y(:, t) .* diag(0.1 * y2) .* diag(y2);
        case 3
            g = y(:, t) .* sqrt(diag(y2));
        case 4
            g = y(:, t) .* log(diag(y2));
        case 5
            g = -y(:, t) .* (diag(y2) < 0.9);
        case 6
            g = y(:, t) ./ log(diag(y2));
        case 7
            g = -y(:, t) ./ sqrt(diag(y2));
        case 8
            g = -y(:, t) ./ diag(y2);
    end
    
    gy = g * y(:, t)'; % Inner product of the nonlinearity adaptation and the separated source
    G = (y2 - idsou) / (1 + lambda * trace(y2)) + (gy - gy') / (1 + lambda * abs(g' * y(:, t))); % Update matrix G
    B = B - lambda * G * B; % Update the separation matrix using matrix G
    conv(t) = norm(G); % Store the convergence value at each iteration
end
