function [penalty_factor cc ee L_curveC L_curveE] = OptimalSmoothnessFactor(xx, bb, ss, mode, SmoothnessFactor, penalty_factor0, signorm, ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator)
if(mode == 0)
    penalty_factor = SmoothnessFactor; % lambda
    cc = Cost_lambda(xx, bb, ss, 0, penalty_factor);
    ee = Cost_gamma(xx, bb, ss, 0, 1/penalty_factor);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
elseif(mode == 1)
    penalty_factor = 1/SmoothnessFactor; % gamma
    cc = Cost_gamma(xx, bb, ss, 0, penalty_factor)/length(xx);
    ee = Cost_lambda(xx, bb, ss, 0, 1/penalty_factor);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
elseif(mode == 2)
    if(Cost_lambda(xx, bb, ss, SmoothnessFactor, 0) <= 0)
        lambda = 0;
    else
        min_lambda = 0;
        if(isinf(penalty_factor0))
            max_lambda = 1;
        elseif(penalty_factor0 >= ACCURACY)
            max_lambda = penalty_factor0;
        else
            max_lambda = 1;
        end
        while(Cost_lambda(xx, bb, ss, SmoothnessFactor, max_lambda) > 0)
            min_lambda = max_lambda;
            max_lambda = 10*max_lambda;
        end
        % bisection
        iteration_counter = 0;
        while(abs(max_lambda - min_lambda) > ACCURACY  && iteration_counter < MAX_ITERATION)
            mid_lambda = (min_lambda + max_lambda)/2;
            c_min = Cost_lambda(xx, bb, ss, SmoothnessFactor, min_lambda);
            c_mid = Cost_lambda(xx, bb, ss, SmoothnessFactor, mid_lambda);
            if(c_min * c_mid <= 0)
                max_lambda = mid_lambda;
            else
                min_lambda = mid_lambda;
            end
            iteration_counter = iteration_counter + 1;
        end
        lambda = (min_lambda + max_lambda)/2;
    end
    penalty_factor = lambda;
    if(~isinf(penalty_factor0))
        penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
    end
    cc = Cost_lambda(xx, bb, ss, 0, lambda);
    ee = Cost_gamma(xx, bb, ss, 0, 1/lambda);
    sweep = linspace(lambda/LCurveSweepFator, LCurveSweepFator*lambda, LCurveSweepLength);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
    for pp = 1 : LCurveSweepLength,
        L_curveC(pp) = Cost_lambda(xx, bb, ss, 0, sweep(pp));
        L_curveE(pp) = Cost_gamma(xx, bb, ss, 0, 1/sweep(pp));
    end
elseif(mode == 3)
    if(Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, 0) <= 0)
        gamma = 0;
    else
        min_gamma = 0;
        if(isinf(penalty_factor0))
            max_gamma = 1;
        elseif(penalty_factor0 >= ACCURACY)
            max_gamma = penalty_factor0;
        else
            max_gamma = 1;
        end
        while(Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, max_gamma) > 0)
            min_gamma = max_gamma;
            max_gamma = 10*max_gamma;
        end
        % bisection
        iteration_counter = 0;
        while(abs(max_gamma - min_gamma) > ACCURACY && iteration_counter < MAX_ITERATION)
            mid_gamma = (min_gamma + max_gamma)/2;
            c_min = Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, min_gamma);
            c_mid = Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, mid_gamma);
            if(c_min * c_mid <= 0)
                max_gamma = mid_gamma;
            else
                min_gamma = mid_gamma;
            end
            iteration_counter = iteration_counter + 1;
        end
        gamma = (min_gamma + max_gamma)/2;
    end
    penalty_factor = 1/gamma; % lambda
    if(~isinf(penalty_factor0))
        penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
    end
    cc = Cost_gamma(xx, bb, ss, 0, gamma)/length(xx);
    ee = Cost_lambda(xx, bb, ss, 0, 1/gamma);
    sweep = linspace(gamma/LCurveSweepFator, LCurveSweepFator*gamma, LCurveSweepLength);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
    for pp = 1 : LCurveSweepLength,
        L_curveC(pp) = Cost_gamma(xx, bb, ss, 0, sweep(pp))/length(xx);
        L_curveE(pp) = Cost_lambda(xx, bb, ss, 0, 1/sweep(pp));
    end
elseif(mode == 4) % naive lambda adaptation
    penalty_factor = signorm*sqrt(length(xx))/norm(xx)*SmoothnessFactor;%(gamma0/(sqrt(N/Nk)*norm(xx)/norm(x)));
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 5) % naive gamma adaptation
    penalty_factor = 1/signorm*sqrt(1/length(xx))*norm(xx)*SmoothnessFactor;%(gamma0/(sqrt(N/Nk)*norm(xx)/norm(x)));
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 6) % L-curve based lambda
    if(Cost_lambda(xx, bb, ss, SmoothnessFactor, 0) <= 0)
        lambda = 0;
    else
        min_lambda = 0;
        if(isinf(penalty_factor0))
            max_lambda = 1;
        elseif(penalty_factor0 >= ACCURACY)
            max_lambda = penalty_factor0;
        else
            max_lambda = 1;
        end
        while(Cost_lambda(xx, bb, ss, SmoothnessFactor, max_lambda) > 0)
            min_lambda = max_lambda;
            max_lambda = 10*max_lambda;
        end
        % bisection
        iteration_counter = 0;
        while(abs(max_lambda - min_lambda) > ACCURACY  && iteration_counter < MAX_ITERATION)
            mid_lambda = (min_lambda + max_lambda)/2;
            c_min = Cost_lambda(xx, bb, ss, SmoothnessFactor, min_lambda);
            c_mid = Cost_lambda(xx, bb, ss, SmoothnessFactor, mid_lambda);
            if(c_min * c_mid <= 0)
                max_lambda = mid_lambda;
            else
                min_lambda = mid_lambda;
            end
            iteration_counter = iteration_counter + 1;
        end
        lambda = (min_lambda + max_lambda)/2;
    end
    cc = Cost_lambda(xx, bb, ss, 0, lambda);
    ee = Cost_gamma(xx, bb, ss, 0, 1/lambda);
    sweep = linspace(lambda/LCurveSweepFator, LCurveSweepFator*lambda, LCurveSweepLength);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
    for pp = 1 : LCurveSweepLength,
        L_curveC(pp) = Cost_lambda(xx, bb, ss, 0, sweep(pp));
        L_curveE(pp) = Cost_gamma(xx, bb, ss, 0, 1/sweep(pp));
    end
    sm = L_curveC + L_curveC;
    [~, ind] = min(sm);
    penalty_factor = sweep(ind);
    if(~isinf(penalty_factor0))
        penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
    end
elseif(mode == 7) % L-curve based gamma
    if(Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, 0) <= 0)
        gamma = 0;
    else
        min_gamma = 0;
        if(isinf(penalty_factor0))
            max_gamma = 1;
        elseif(penalty_factor0 >= ACCURACY)
            max_gamma = penalty_factor0;
        else
            max_gamma = 1;
        end
        while(Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, max_gamma) > 0)
            min_gamma = max_gamma;
            max_gamma = 10*max_gamma;
        end
        % bisection
        iteration_counter = 0;
        while(abs(max_gamma - min_gamma) > ACCURACY && iteration_counter < MAX_ITERATION)
            mid_gamma = (min_gamma + max_gamma)/2;
            c_min = Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, min_gamma);
            c_mid = Cost_gamma(xx, bb, ss, length(xx)*SmoothnessFactor, mid_gamma);
            if(c_min * c_mid <= 0)
                max_gamma = mid_gamma;
            else
                min_gamma = mid_gamma;
            end
            iteration_counter = iteration_counter + 1;
        end
        gamma = (min_gamma + max_gamma)/2;
    end
    cc = Cost_gamma(xx, bb, ss, 0, gamma)/length(xx);
    ee = Cost_lambda(xx, bb, ss, 0, 1/gamma);
    sweep = linspace(gamma/LCurveSweepFator, LCurveSweepFator*gamma, LCurveSweepLength);
    L_curveC = zeros(1, LCurveSweepLength);
    L_curveE = zeros(1, LCurveSweepLength);
    for pp = 1 : LCurveSweepLength,
        L_curveC(pp) = Cost_gamma(xx, bb, ss, 0, sweep(pp))/length(xx);
        L_curveE(pp) = Cost_lambda(xx, bb, ss, 0, 1/sweep(pp));
    end
    sm = L_curveC + L_curveC;
    [~, ind] = min(sm);
    penalty_factor = 1/sweep(ind); % lambda
    if(~isinf(penalty_factor0))
        penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
    end
else
    penalty_factor = SmoothnessFactor;
    %     warning('Invalid mode; using default mode = 0');
end
