function [x_smoothed1, x_smoothed2, x_smoothed3, optim_gammas1, optim_gammas2] = ECGSmoothnessPriorsDenoiserBWFiltFilt(x, SmoothnessFactor, varargin)
% Block-wise ECG denoising based on smoothness priors
%
% [x_smoothed1, x_smoothed2] = ECGSmoothnessPriorsDenoiserBW(x, SmoothnessFactor, mode, FilterParam, KnotsParam, ACCURACY, MAX_ITERATION)
%
% mode = 0: y_opt = argmin(|x - y| + gamma*|D*y + b|), fixed smoothness penalty (gamma)
% mode = 1: y_opt = argmin(|D * y + b| + lambda*|x - y|), fixed MSE penalty (lambda)
% mode = 2: y_opt = argmin(|x - y|), s.t. |D*y + b| = epsilon, fixed smoothness (epsilon)
% mode = 3: y_opt = argmin(|D * y + b|), s.t. |x - y| = N*sigma, fixed noise variance (sigma)
%
% Note: Segment lengths should not be shorter than the filter length
% Reza Sameni, 2015
%

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


% FILTFILT VERSION
phi = conv(h, h(end:-1:1));
% % % phi = conv(h, h);
gg = phi*SmoothnessFactor;
midindex = (length(phi)+1)/2;
gg(midindex) = gg(midindex) + 1;
r = roots(gg);
r_abs = abs(r);
I_causal = r_abs < 1;%0.9999;
r_causal = r(I_causal);
p_causal = poly(r_causal);
% I_anticausal = r_abs > 1;%0.9999;
% r_anticausal = r(I_anticausal);
% p_anticausal = poly(r_anticausal);

M = size(x,1);
N = size(x,2);
L = length(h);
guardlen_l = L - 1; % left guard window lengh, (L - 1) is enough!
guardlen_r = L - 1; % right guard window lengh, (L - 1) is enough!
% orvlp = 1; % segment overalp should be 0; but one sample of overlap is added to the end of each segment to improve the continuity for the first run
orvlp = 0; % segment overalp should be 0; but one sample of overlap is added to the end of each segment to improve the continuity for the first run
x_smoothed1 = zeros(size(x));
x_smoothed2 = zeros(size(x));
x_smoothed3 = zeros(size(x));

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
gamma = 0;
Dd_pre = [];
% optim_gammas1 = zeros(M, length(knots0));
optim_gammas1 = zeros(size(x));
for k = 1 : length(knots0) - 1,
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
    
    for ch = 1 : M,
        % take a segment
        xx = x(ch, indexes)';
        
        if(~isequal(Dd, Dd_pre))
            [~, S, V] = svd(Dd);
            Dd_pre = Dd;
        end
        xx_ = V.'*xx;
        b_ = zeros(size(Dd, 1), 1);
        s_ = diag(S);
        % regularization factor update
        gamma = OptimalSmoothnessFactor(xx_, b_, s_, mode, SmoothnessFactor, gamma, norm(x)/sqrt(N), ACCURACY, MAX_ITERATION, FORGETTING_FACTOR);
        %         optim_gammas1(ch, k) = gamma;
        optim_gammas1(ch, indexes) = gamma;
        
        % smoothing
        if(isinf(gamma))
            x_smoothed1(ch, indexes) = xx;
            gamma = 0;
        else
            x_smoothed1(ch, indexes) = (eye(size(Dd, 2)) + gamma*(Dd'*Dd))\xx;
        end
    end
end

D_pre = [];
% optim_gammas2 = zeros(M, length(knots));
optim_gammas2 = zeros(size(x));
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
        
        % boundary vectors from the first round
        l = x_smoothed1(ch, indxl)';
        r = x_smoothed1(ch, indxr)';
        b = D_l*l + D_r*r;
        
        if(~isequal(D, D_pre))
            [U, S, V] = svd(D);
            D_pre = D;
        end
        xx_ = V.'*xx;
        b_ = U.'*b;
        s_ = diag(S);
        % regularization factor update
        gamma = OptimalSmoothnessFactor(xx_, b_, s_, mode, SmoothnessFactor, gamma, norm(x)/sqrt(N), ACCURACY, MAX_ITERATION, FORGETTING_FACTOR);
        %         optim_gammas2(ch, k) = gamma;
        optim_gammas2(ch, indexes) = gamma;
        
        % smoothing
        if(isinf(gamma))
            x_smoothed2(ch, indexes) = xx;
            x_smoothed3(ch, indexes) = xx;
            gamma = 0;
        else
            x_smoothed2(ch, indexes) = (eye(size(D, 2)) + gamma*(D'*D))\(xx - gamma*D'*b);
%             forward = filter(sum(p_causal), p_causal, xx, zi);
            %             backward = filter(sum(p_causal), p_causal, forward(end : -1 : 1), r(end : -1 : 1));
            %             backward = filter(sum(p_anticausal), p_anticausal, forward(end : -1 : 1), r(end : -1 : 1));
            forward = [l ; zeros(length(indexes),1)];
            for m = 1 : length(indexes),
                for p = 2 : length(p_causal),
                    forward(m + length(l)) = forward(m + length(l)) - p_causal(p)*forward(m + length(l) - (p-1));
                end
                forward(m + length(l)) = (forward(m + length(l)) + sum(p_causal)*xx(m)) / p_causal(1);
            end
            forward = forward(length(l)+1 : end);
            
            % flip the signal
            forward = forward(end: -1 : 1);
            
            backward = [r(end:-1:1) ; zeros(length(indexes),1)];
            for m = 1 : length(indexes),
                for p = 2 : length(p_causal),
                    backward(m + length(r)) = backward(m + length(r)) - p_causal(p)*backward(m + length(r) - (p-1));
                end
                backward(m + length(r)) = (backward(m + length(r)) + sum(p_causal)*forward(m)) / p_causal(1);
            end
            backward = backward(length(r)+1 : end);
            
            % flip back the signal
            backward = backward(end: -1 : 1);
            
            x_smoothed3(ch, indexes) = backward;
        end
    end
end