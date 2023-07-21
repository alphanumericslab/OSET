function [penalty_factor, cc, ee, L_curveC, L_curveE] = optimal_smoothness_factor(xx, bb, ss, mode, SmoothnessFactor, penalty_factor0, signorm, ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator)
% optimal_smoothness_factor - Calculate the optimal smoothness factor for a piece-wise Tikhonov regularization algorithm.
%
% Usage:
%   [penalty_factor, cc, ee, L_curveC, L_curveE] = optimal_smoothness_factor(xx, bb, ss, mode, SmoothnessFactor, penalty_factor0, signorm, ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator)
%
% Inputs:
%   xx: Input signal.
%   bb: Signal prior.
%   ss: Smoothness prior.
%   mode: Mode of operation (0 to 7).
%       - Mode 0: Direct Smoothness Factor (lambda)
%       - Mode 1: Inverse Smoothness Factor (gamma)
%       - Mode 2: Bisection Method for Optimal Smoothness Factor (lambda)
%       - Mode 3: Bisection Method for Optimal Smoothness Factor (gamma)
%       - Mode 4: Naive lambda adaptation
%       - Mode 5: Naive gamma adaptation
%       - Mode 6: L-curve based lambda
%       - Mode 7: L-curve based gamma
%   SmoothnessFactor: The user-specified smoothness factor.
%   penalty_factor0: Initial penalty factor value.
%   signorm: Scale factor for the penalty factor calculation.
%   ACCURACY: Convergence accuracy for the bisection method.
%   MAX_ITERATION: Maximum number of iterations for the bisection method.
%   FORGETTING_FACTOR: Forgetting factor for the penalty factor.
%   LCurveSweepLength: Length of the L-curve sweep.
%   LCurveSweepFator: L-curve sweep factor.
%
% Outputs:
%   penalty_factor: Optimal smoothness factor (penalty factor) value.
%   cc: Cost value for penalty_factor.
%   ee: Cost value for 1/penalty_factor.
%   L_curveC: L-curve values for penalty_factor.
%   L_curveE: L-curve values for 1/penalty_factor.
%
% Reference: Sameni, R. (2017). Online filtering using piecewise smoothness
%   priors: Application to normal and abnormal electrocardiogram denoising.
%   In Signal Processing (Vol. 133, pp. 52–63).
%   https://doi.org/10.1016/j.sigpro.2016.10.019
%
% Revision History:
%   2015: First release
%   2023: Documented and renamed from deprecated version OptimalSmoothnessFactor
%
% Reza Sameni, 2015-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET


switch mode
    case 0
        penalty_factor = SmoothnessFactor; % lambda
        cc = cost_lambda(xx, bb, ss, 0, penalty_factor);
        ee = cost_gamma(xx, bb, ss, 0, 1/penalty_factor);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
    case 1
        penalty_factor = 1/SmoothnessFactor; % gamma
        cc = cost_gamma(xx, bb, ss, 0, penalty_factor)/length(xx);
        ee = cost_lambda(xx, bb, ss, 0, 1/penalty_factor);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
    case 2
        if (cost_lambda(xx, bb, ss, SmoothnessFactor, 0) <= 0)
            lambda = 0;
        else
            min_lambda = 0;
            if (isinf(penalty_factor0))
                max_lambda = 1;
            elseif (penalty_factor0 >= ACCURACY)
                max_lambda = penalty_factor0;
            else
                max_lambda = 1;
            end
            while (cost_lambda(xx, bb, ss, SmoothnessFactor, max_lambda) > 0)
                min_lambda = max_lambda;
                max_lambda = 10 * max_lambda;
            end
            % bisection
            iteration_counter = 0;
            while (abs(max_lambda - min_lambda) > ACCURACY  && iteration_counter < MAX_ITERATION)
                mid_lambda = (min_lambda + max_lambda) / 2;
                c_min = cost_lambda(xx, bb, ss, SmoothnessFactor, min_lambda);
                c_mid = cost_lambda(xx, bb, ss, SmoothnessFactor, mid_lambda);
                if (c_min * c_mid <= 0)
                    max_lambda = mid_lambda;
                else
                    min_lambda = mid_lambda;
                end
                iteration_counter = iteration_counter + 1;
            end
            lambda = (min_lambda + max_lambda) / 2;
        end
        penalty_factor = lambda;
        if (~isinf(penalty_factor0))
            penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
        end
        cc = cost_lambda(xx, bb, ss, 0, lambda);
        ee = cost_gamma(xx, bb, ss, 0, 1/lambda);
        sweep = linspace(lambda / LCurveSweepFator, LCurveSweepFator * lambda, LCurveSweepLength);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
        for pp = 1 : LCurveSweepLength
            L_curveC(pp) = cost_lambda(xx, bb, ss, 0, sweep(pp));
            L_curveE(pp) = cost_gamma(xx, bb, ss, 0, 1 / sweep(pp));
        end
    case 3
        if (cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, 0) <= 0)
            gamma = 0;
        else
            min_gamma = 0;
            if (isinf(penalty_factor0))
                max_gamma = 1;
            elseif (penalty_factor0 >= ACCURACY)
                max_gamma = penalty_factor0;
            else
                max_gamma = 1;
            end
            while (cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, max_gamma) > 0)
                min_gamma = max_gamma;
                max_gamma = 10 * max_gamma;
            end
            % bisection
            iteration_counter = 0;
            while (abs(max_gamma - min_gamma) > ACCURACY && iteration_counter < MAX_ITERATION)
                mid_gamma = (min_gamma + max_gamma) / 2;
                c_min = cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, min_gamma);
                c_mid = cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, mid_gamma);
                if (c_min * c_mid <= 0)
                    max_gamma = mid_gamma;
                else
                    min_gamma = mid_gamma;
                end
                iteration_counter = iteration_counter + 1;
            end
            gamma = (min_gamma + max_gamma) / 2;
        end
        penalty_factor = 1 / gamma; % lambda
        if (~isinf(penalty_factor0))
            penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
        end
        cc = cost_gamma(xx, bb, ss, 0, gamma) / length(xx);
        ee = cost_lambda(xx, bb, ss, 0, 1 / gamma);
        sweep = linspace(gamma / LCurveSweepFator, LCurveSweepFator * gamma, LCurveSweepLength);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
        for pp = 1 : LCurveSweepLength
            L_curveC(pp) = cost_gamma(xx, bb, ss, 0, sweep(pp)) / length(xx);
            L_curveE(pp) = cost_lambda(xx, bb, ss, 0, 1 / sweep(pp));
        end
    case 4 % naive lambda adaptation
        penalty_factor = signorm * sqrt(length(xx)) / norm(xx) * SmoothnessFactor; % (gamma0/(sqrt(N/Nk)*norm(xx)/norm(x)));
        penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
    case 5 % naive gamma adaptation
        penalty_factor = 1 / signorm * sqrt(1 / length(xx)) * norm(xx) * SmoothnessFactor; % (gamma0/(sqrt(N/Nk)*norm(xx)/norm(x)));
        penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
    case 6 % L-curve based lambda
        if (cost_lambda(xx, bb, ss, SmoothnessFactor, 0) <= 0)
            lambda = 0;
        else
            min_lambda = 0;
            if (isinf(penalty_factor0))
                max_lambda = 1;
            elseif (penalty_factor0 >= ACCURACY)
                max_lambda = penalty_factor0;
            else
                max_lambda = 1;
            end
            while (cost_lambda(xx, bb, ss, SmoothnessFactor, max_lambda) > 0)
                min_lambda = max_lambda;
                max_lambda = 10 * max_lambda;
            end
            % bisection
            iteration_counter = 0;
            while (abs(max_lambda - min_lambda) > ACCURACY  && iteration_counter < MAX_ITERATION)
                mid_lambda = (min_lambda + max_lambda) / 2;
                c_min = cost_lambda(xx, bb, ss, SmoothnessFactor, min_lambda);
                c_mid = cost_lambda(xx, bb, ss, SmoothnessFactor, mid_lambda);
                if (c_min * c_mid <= 0)
                    max_lambda = mid_lambda;
                else
                    min_lambda = mid_lambda;
                end
                iteration_counter = iteration_counter + 1;
            end
            lambda = (min_lambda + max_lambda) / 2;
        end
        cc = cost_lambda(xx, bb, ss, 0, lambda);
        ee = cost_gamma(xx, bb, ss, 0, 1 / lambda);
        sweep = linspace(lambda / LCurveSweepFator, LCurveSweepFator * lambda, LCurveSweepLength);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
        for pp = 1 : LCurveSweepLength
            L_curveC(pp) = cost_lambda(xx, bb, ss, 0, sweep(pp));
            L_curveE(pp) = cost_gamma(xx, bb, ss, 0, 1 / sweep(pp));
        end
        sm = L_curveC + L_curveC;
        [~, ind] = min(sm);
        penalty_factor = sweep(ind);
        if (~isinf(penalty_factor0))
            penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
        end
    case 7 % L-curve based gamma
        if (cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, 0) <= 0)
            gamma = 0;
        else
            min_gamma = 0;
            if (isinf(penalty_factor0))
                max_gamma = 1;
            elseif (penalty_factor0 >= ACCURACY)
                max_gamma = penalty_factor0;
            else
                max_gamma = 1;
            end
            while (cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, max_gamma) > 0)
                min_gamma = max_gamma;
                max_gamma = 10 * max_gamma;
            end
            % bisection
            iteration_counter = 0;
            while (abs(max_gamma - min_gamma) > ACCURACY && iteration_counter < MAX_ITERATION)
                mid_gamma = (min_gamma + max_gamma) / 2;
                c_min = cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, min_gamma);
                c_mid = cost_gamma(xx, bb, ss, length(xx) * SmoothnessFactor, mid_gamma);
                if (c_min * c_mid <= 0)
                    max_gamma = mid_gamma;
                else
                    min_gamma = mid_gamma;
                end
                iteration_counter = iteration_counter + 1;
            end
            gamma = (min_gamma + max_gamma) / 2;
        end
        penalty_factor = 1 / gamma; % lambda
        if (~isinf(penalty_factor0))
            penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
        end
        cc = cost_gamma(xx, bb, ss, 0, gamma) / length(xx);
        ee = cost_lambda(xx, bb, ss, 0, 1 / gamma);
        sweep = linspace(gamma / LCurveSweepFator, LCurveSweepFator * gamma, LCurveSweepLength);
        L_curveC = zeros(1, LCurveSweepLength);
        L_curveE = zeros(1, LCurveSweepLength);
        for pp = 1 : LCurveSweepLength
            L_curveC(pp) = cost_gamma(xx, bb, ss, 0, sweep(pp)) / length(xx);
            L_curveE(pp) = cost_lambda(xx, bb, ss, 0, 1 / sweep(pp));
        end
        sm = L_curveC + L_curveC;
        [~, ind] = min(sm);
        penalty_factor = 1 / sweep(ind); % lambda
        if (~isinf(penalty_factor0))
            penalty_factor = FORGETTING_FACTOR * penalty_factor + (1 - FORGETTING_FACTOR) * penalty_factor0;
        end
    otherwise
        error('Invalid mode');
end

end

%//////////////////////////////////////////////////////////////////////////
function C = cost_gamma(x, b, s, sigma2, gamma)
n = min([length(x), length(b), length(s)]);
x = x(1:n);
b = b(1:n);
s = s(1:n);
C = sum(((s.^2.*x + s.*b)./(gamma + s.^2)).^2) - sigma2;
end

%//////////////////////////////////////////////////////////////////////////
function C = cost_lambda(x, b, s, epsilon2, lambda)
n = min([length(x), length(b), length(s)]);
x = x(1:n);
b = b(1:n);
s = s(1:n);
C = sum(((s.*x + b)./(1 + lambda.*s.^2)).^2) - epsilon2;
end
