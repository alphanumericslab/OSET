function [x_smoothed1, x_smoothed2, optim_lambdas1, optim_lambdas2, cc1, ee1, cc2, ee2, L_curveC1, L_curveE1, L_curveC2, L_curveE2, knots0, knots] = ECGSmoothnessPriorsDenoiserBW(x, SmoothnessFactor, varargin)
% Block-wise ECG denoising based on smoothness priors
%
% [x_smoothed1, x_smoothed2] = ECGSmoothnessPriorsDenoiserBW(x, SmoothnessFactor, mode, FilterParam, KnotsParam, ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength)
%
% mode = 0: y_opt = argmin(|x - y| + lambda*|D*y + b|), fixed smoothness penalty (lambda)
% mode = 1: y_opt = argmin(|D * y + b| + gamma*|x - y|), fixed MSE penalty (gamma)
% mode = 2: y_opt = argmin(|x - y|), s.t. |D*y + b| = epsilon, fixed smoothness (epsilon)
% mode = 3: y_opt = argmin(|D * y + b|), s.t. |x - y| = N*sigma, fixed noise variance (sigma)
% mode = 4: y_opt = argmin(|x - y| + lambda*|D*y + b|), fixed smoothness penalty (naive lambda adaptation)
% mode = 5: y_opt = argmin(|D * y + b| + gamma*|x - y|), fixed MSE penalty (naive gamma adaptation)
% mode = 6: y_opt = argmin(|x - y| + lambda*|D*y + b|), fixed smoothness penalty (L-curve based lambda estimation)
% mode = 7: y_opt = argmin(|D * y + b| + gamma*|x - y|), fixed MSE penalty (L-curve based lambda gamma estimation)
%
% Note: Segment lengths should not be shorter than the filter length
%
% Reza Sameni, Copyright 2015
% email: reza.sameni@gmail.com

if(nargin > 2 && ~isempty(varargin{1}))
    mode = varargin{1};
else
    mode = 0;
end

if(nargin > 3)
    FilterParam = varargin{2};
else
    FilterParam = [];
end

if(nargin > 4)
    KnotsParam = varargin{3};
else
    KnotsParam = [];
end

if(nargin > 5 && ~isempty(varargin{4}))
    ACCURACY = varargin{4};
else
    ACCURACY = 1e-8;
end

if(nargin > 6 && ~isempty(varargin{5}))
    MAX_ITERATION = varargin{5};
else
    MAX_ITERATION = 250;
end

if(nargin > 7 && ~isempty(varargin{6}))
    FORGETTING_FACTOR = varargin{6};
else
    FORGETTING_FACTOR = 1.0; % forget the previous blocks
end

if(nargin > 8 && ~isempty(varargin{7}))
    LCurveSweepLength = varargin{7}; % 500 is good as default
else
    LCurveSweepLength = 0; % No L-curves
end

if(nargin > 9 && ~isempty(varargin{8}))
    LCurveSweepFator = varargin{8}; % 500 is good as default
else
    LCurveSweepFator = 25;
end


% use the given filter or calculate the filter impulse response
if(isempty(FilterParam))
    DiffOrder = 2;
    h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
elseif(length(FilterParam) == 1)
    DiffOrder = FilterParam;
    h = diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder);
else
    h = FilterParam(:)';
end

M = size(x,1);
N = size(x,2);
L = length(h);
guardlen_l = L - 1; % left guard window lengh, (L - 1) is enough!
guardlen_r = L - 1; % right guard window lengh, (L - 1) is enough!
orvlp = 0; % segment overalp should be 0
x_smoothed1 = zeros(size(x));
x_smoothed2 = zeros(size(x));

% use the user-defined knots or calculate the knot points
if(isempty(KnotsParam))
    wlen = round(N/10); % 10 segments by default, if not specified by user
    knots = 1 : round(wlen) : N + 1;
    D0 = toeplitz([h(1) zeros(1, wlen - L + orvlp)], [h zeros(1, wlen - L + orvlp)]);
    D1 = toeplitz([h(1) zeros(1, wlen + guardlen_l + guardlen_r - L + orvlp)], [h zeros(1, wlen + guardlen_l + guardlen_r - L + orvlp)]);
    %     midindex = round(length(h)/2);
    %     D0 = toeplitz([h(midindex:-1:1) zeros(1, wlen - midindex)], [h(midindex:end) zeros(1, wlen - length(h) + midindex - 1)]);
    %     [U0, S0, V0] = svd(D0);
elseif(length(KnotsParam) == 1)
    wlen = KnotsParam;
    knots = 1 : round(wlen) : N + 1;
    D0 = toeplitz([h(1) zeros(1, wlen - L + orvlp)], [h zeros(1, wlen - L + orvlp)]);
    D1 = toeplitz([h(1) zeros(1, wlen + guardlen_l + guardlen_r - L + orvlp)], [h zeros(1, wlen + guardlen_l + guardlen_r - L + orvlp)]);
    %     midindex = round(length(h)/2);
    %     D0 = toeplitz([h(midindex:-1:1) zeros(1, wlen - midindex)], [h(midindex:end) zeros(1, wlen - length(h) + midindex - 1)]);
    %     [U0, S0, V0] = svd(D0);
else
    wlen = 0;
    knots = KnotsParam;
end

% Preprocess the knots sequence
% remove illegal knots
knots(knots < 1) = [];
knots(knots > (N + 1 - orvlp)) = [];
knots = sort(knots); % sort the knots if they are not sorted
% add the first and last samples as knots to avoid end-point discontinuities
if(knots(1) > 1)
    knots = [1, knots];
end
if(knots(end) < (N + 1 - orvlp))
    knots = [knots, (N + 1 - orvlp)];
end

% Knots for the first round
knots0 = round((knots(1:end-1) + knots(2:end))/2);
if(knots0(1) > 1)
    knots0 = [1, knots0];
end
if(knots0(end) < (N + 1 - orvlp))
    knots0 = [knots0, (N + 1 - orvlp)];
end

% First round of smoothing (without any boundary conditions)
lambda = 0;
Dd_pre = [];
% optim_lambdas1 = zeros(M, length(knots0));
optim_lambdas1 = zeros(size(x));
cc1 = zeros(M, length(knots0)-1);
ee1 = zeros(M, length(knots0)-1);
L_curveC1 = zeros(M, length(knots0)-1, LCurveSweepLength);
L_curveE1 = zeros(M, length(knots0)-1, LCurveSweepLength);
for k = 1 : length(knots0) - 1
    indexes = knots0(k) : knots0(k + 1) - 1 + orvlp;
    if(length(indexes) < L)
        pad = L - length(indexes);
        indexes_extended = [indexes indexes(zeros(1, pad) + length(indexes))];
        indexes = indexes_extended;
    end
    Nk = length(indexes);
    if(Nk == wlen + orvlp)
        Dd = D0;
    else
        Dd = toeplitz([h(1) zeros(1, Nk - L)], [h zeros(1, Nk - L)]);
    end
    
    for ch = 1 : M
        % take a segment
        xx = x(ch, indexes)';
        mn = mean(xx);
%         xx = xx - mn; % USING THE MEAN CAUSES DISCONTINUITIES AT THE SEGMENT BOUNDARIES
        
        if(~isequal(Dd, Dd_pre))
            [~, S, V] = svd(Dd);
            Dd_pre = Dd;
        end
        xx_tilde = V.'*xx;
        b_tilde = zeros(size(Dd, 1), 1);
        s_tilde = diag(S);
        % regularization factor update
        
        % ATTENTION: UNDER EVALUATION
        [lambda cc ee L_curveC L_curveE] = optimal_smoothness_factor(xx_tilde, b_tilde, s_tilde, mode, SmoothnessFactor, lambda, norm(x)/sqrt(N), ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator);
        %         [lambda cc ee L_curveC L_curveE] = OptimalSmoothnessFactor(xx_tilde, b_tilde, s_tilde, 1, SmoothnessFactor, lambda, norm(x)/sqrt(N), ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator);
        
        %%% disp(['SmoothnessFactor = ' num2str(SmoothnessFactor) ', ee = ' num2str(ee) ', cc = ' num2str(cc)]); % check to see if the optimal point is correct
        
        %         optim_lambdas1(ch, k) = lambda;
        optim_lambdas1(ch, indexes) = lambda;
        cc1(ch, k) = cc;
        ee1(ch, k) = ee;
        L_curveC1(ch, k, :) = L_curveC;
        L_curveE1(ch, k, :) = L_curveE;
        
        % smoothing
        if(isinf(lambda))   % use gamma-based formulation setting gamma = 0
            x_smoothed1(ch, indexes) = mn;
            %         elseif(cond((eye(size(Dd, 2)) + lambda*(Dd'*Dd))) > 100)   % use gamma-based formulation setting gamma = 0
            %             gamma = 1/lambda;
            %             x_smoothed1(ch, indexes) = (gamma*eye(size(Dd, 2)) + (Dd'*Dd))\(gamma*xx);
        else                % use lambda-based formulation
            x_smoothed1(ch, indexes) = (eye(size(Dd, 2)) + lambda*(Dd'*Dd))\xx;
            
            % equivalent form (refer to the paper):
            %             n = min(size(S));
            %             x_smoothed1(ch, indexes) = V * [(xx_tilde(1:n)./(1 + lambda*s_tilde(1:n).^2)) ; xx_tilde(n+1:end)];
        end
    end
end

D_pre = [];
% optim_lambdas2 = zeros(M, length(knots));
optim_lambdas2 = zeros(size(x));
cc2 = zeros(M, length(knots)-1);
ee2 = zeros(M, length(knots)-1);
L_curveC2 = zeros(M, length(knots)-1, LCurveSweepLength);
L_curveE2 = zeros(M, length(knots)-1, LCurveSweepLength);
% Second round of smoothing (with boundary conditions from the first round)
for k = 1 : length(knots) - 1,
    indexes = knots(k) : knots(k + 1) - 1 + orvlp;
    Nk = length(indexes);
    if(Nk == wlen + orvlp)
        Dd = D1;
    else
        Dd = toeplitz([h(1) zeros(1, Nk + guardlen_l + guardlen_r - L)], [h zeros(1, Nk + guardlen_l + guardlen_r - L)]);
    end
    
    D_l = Dd(:, 1 : guardlen_l);
    D = Dd(:, guardlen_l + 1 : end - guardlen_r);
    D_r = Dd(:, end - guardlen_r + 1 : end);
    
    indxl = indexes(1) - guardlen_l : indexes(1) - 1;
    indxl(indxl < 1) = 1; % prevent out of range indexes
    
    indxr = indexes(end) + 1 : indexes(end) + guardlen_r;
    indxr(indxr > N) = N; % prevent out of range indexes
    
    for ch = 1 : M,
        % take a segment
        xx = x(ch, indexes)';
%         mn = mean(xx);
%         xx = xx - mn;% USING THE MEAN CAUSES DISCONTINUITIES AT THE SEGMENT BOUNDARIES
        
        % boundary vectors from the first round
        l = x_smoothed1(ch, indxl)';
        r = x_smoothed1(ch, indxr)';
        b = D_l*l + D_r*r;
        
        if(~isequal(D, D_pre))
            [U, S, V] = svd(D);
            D_pre = D;
        end
        xx_tilde = V.'*xx;
        b_tilde = U.'*b;
        s_tilde = diag(S);
        % regularization factor update
        [lambda cc ee L_curveC L_curveE] = optimal_smoothness_factor(xx_tilde, b_tilde, s_tilde, mode, SmoothnessFactor, lambda, norm(x)/sqrt(N), ACCURACY, MAX_ITERATION, FORGETTING_FACTOR, LCurveSweepLength, LCurveSweepFator);
        %         disp(['SmoothnessFactor = ' num2str(SmoothnessFactor) ', ee = ' num2str(ee) ', cc = ' num2str(cc)]); % check to see if the optimal point is correct
        
        %         optim_lambdas2(ch, k) = lambda;
        optim_lambdas2(ch, indexes) = lambda;
        cc2(ch, k) = cc;
        ee2(ch, k) = ee;
        L_curveC2(ch, k, :) = L_curveC;
        L_curveE2(ch, k, :) = L_curveE;
        % smoothing
        if(isinf(lambda))   % use gamma-based formulation setting gamma = 0
            %             x_smoothed2(ch, indexes) = -(D'*D)\(D'*b);
            x_smoothed2(ch, indexes) = -D\b;
            %             disp(['Hi ' num2str(indexes(1))]);
            %         elseif(cond((eye(size(D, 2)) + lambda*(D'*D))) > 100)
            %             gamma = 1/lambda;
            %             x_smoothed2(ch, indexes) = (gamma*eye(size(D, 2)) + (D'*D))\(gamma*xx - D'*b);
        else                % use lambda-based formulation
            x_smoothed2(ch, indexes) = (eye(size(D, 2)) + lambda*(D'*D))\(xx - lambda*D'*b);
            
            % equivalent form (refer to the paper):
            %             n = min(size(S));
            %             x_smoothed2(ch, indexes) = V * ((xx_tilde(1:n) - lambda*s_tilde(1:n).*b_tilde(1:n))./(1 + lambda*s_tilde(1:n).^2));
        end
    end
end