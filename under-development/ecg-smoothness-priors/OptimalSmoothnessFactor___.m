function penalty_factor = OptimalSmoothnessFactor(xx, bb, ss, mode, SmoothnessFactor, penalty_factor0, signorm, ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, penalty_factor00)
if(mode == 0)
    penalty_factor = SmoothnessFactor; % gamma
elseif(mode == 1)
    penalty_factor = 1/SmoothnessFactor; % lambda
elseif(mode == 2)
    if(Cost_gamma(xx, bb, ss, SmoothnessFactor, 0) <= 0)
        gamma = 0;
    else
        min_gamma = 0;
        if(penalty_factor0 >= ACCURACY)
            max_gamma = penalty_factor0;
        else
            max_gamma = 1;
        end
        while(Cost_gamma(xx, bb, ss, SmoothnessFactor, max_gamma) > 0)
            min_gamma = max_gamma;
            max_gamma = 10*max_gamma;
        end
        % bisection
        iteration_counter = 0;
        while(abs(max_gamma - min_gamma) > ACCURACY  && iteration_counter < MAX_ITERATION)
            mid_gamma = (min_gamma + max_gamma)/2;
            c_min = Cost_gamma(xx, bb, ss, SmoothnessFactor, min_gamma);
            c_mid = Cost_gamma(xx, bb, ss, SmoothnessFactor, mid_gamma);
            if(c_min * c_mid <= 0)
                max_gamma = mid_gamma;
            else
                min_gamma = mid_gamma;
            end
            iteration_counter = iteration_counter + 1;
        end
        gamma = (min_gamma + max_gamma)/2;
    end
    penalty_factor = gamma;
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 3)
    if(Cost_lambda(xx, bb, ss, SmoothnessFactor, 0) <= 0)
        lambda = 0;
    else
        min_lambda = 0;
        if(penalty_factor0 >= ACCURACY)
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
        while(abs(max_lambda - min_lambda) > ACCURACY && iteration_counter < MAX_ITERATION)
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
    penalty_factor = 1/lambda; % gamma
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 4) % naive gamma adaptation
    penalty_factor = signorm*sqrt(length(xx))/norm(xx)*SmoothnessFactor;%(lambda0/(sqrt(N/Nk)*norm(xx)/norm(x)));
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 5) % naive lambda adaptation
    penalty_factor = 1/signorm*sqrt(1/length(xx))*norm(xx)*SmoothnessFactor;%(lambda0/(sqrt(N/Nk)*norm(xx)/norm(x)));
    penalty_factor = FORGETTING_FACTOR*penalty_factor + (1 - FORGETTING_FACTOR)*penalty_factor0;
elseif(mode == 6) % adaptive gamma adaptation
    e0 = Cost_gamma(xx, bb, ss, SmoothnessFactor, penalty_factor0);
    e00 = Cost_gamma(xx, bb, ss, SmoothnessFactor, penalty_factor00);
    if(abs(e0) < ACCURACY)
        penalty_factor = penalty_factor0;
    elseif(abs(e00) < ACCURACY)
        penalty_factor = penalty_factor00;
    else
        q = e0/e00;
        penalty_factor = 1/(1+q)*penalty_factor0 + q/(1+q)*penalty_factor00;
    end
elseif(mode == 7) % adaptive lambda adaptation
    lambda0 = 1/penalty_factor0;
    lambda00 = 1/penalty_factor00;
    %     e0 = Cost_lambda(xx, bb, ss, SmoothnessFactor, lambda0);
    %     e00 = Cost_lambda(xx, bb, ss, SmoothnessFactor, lambda00);
    %     if(abs(e0) < ACCURACY)
    %         penalty_factor = penalty_factor0;
    %     elseif(abs(e00) < ACCURACY)
    %         penalty_factor = penalty_factor00;
    %     else
    %         q = e0/e00;
    %         lambda = 1/(1+q)*lambda0 + q/(1+q)*lambda00;
    %         penalty_factor = 1/lambda; % gamma
    %     end
    c0 = Cost_lambda(xx, bb, ss, 0, lambda0);
    c00 = Cost_lambda(xx, bb, ss, 0, lambda00);
    alpha = (SmoothnessFactor  - c00)/(c0 - c00);
    lambda = alpha*lambda0 + (1 - alpha)*lambda00;
    lambda = max(0, lambda);
    penalty_factor = 1/lambda; % gamma
else
    penalty_factor = SmoothnessFactor;
    %     warning('Invalid mode; using default mode = 0');
end
end

%//////////////////////////////////////////////////////////////////////////
function C = Cost_gamma(x, b, s, epsilon2, gamma)
n = min([length(x), length(b), length(s)]);
x = x(1:n);
% b = x(1:n);
b = b(1:n);
s = s(1:n);
C = sum(((s.*x + b)./(1 + gamma.*s.^2)).^2) - epsilon2;
end

%//////////////////////////////////////////////////////////////////////////
function C = Cost_lambda(x, b, s, sigma2, lambda)
n = min([length(x), length(b), length(s)]);
x = x(1:n);
% b = x(1:n);
b = b(1:n);
s = s(1:n);
C = mean(((s.*b + s.^2.*x)./(lambda + s.^2)).^2) - sigma2;
end
