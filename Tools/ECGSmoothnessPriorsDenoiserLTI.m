function y = ECGSmoothnessPriorsDenoiserLTI(x, SmoothnessFactor, mode, FilterParam)
% LTI ECG denoising based on smoothness priors
%
% [x_smoothed1, x_smoothed2] = ECGSmoothnessPriorsDenoiserBW(x, FilterParam, KnotsParam, mode, lambda0, adapt)
%
% mode = 0: y_opt = argmin(|x - y| + lambda*|D*y + b|), fixed smoothness penalty (lambda)
% mode = 1: y_opt = argmin(|D * y + b| + gamma*|x - y|), fixed MSE penalty (gamma)
%
% Note: Segment lengths should not be shorter than the filter length
% Reza Sameni, 2015
% 
% Modified Oct 2020 to work with multichannel signals

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

th = 0.9999;

phi = conv(h, h(end:-1:1));
    midindex = (length(phi)+1)/2;
if(mode == 0) % lambda-based
    gg = phi*SmoothnessFactor;
    gg(midindex) = gg(midindex) + 1;
    r = roots(gg);
    r_abs = abs(r);
    I_causal = r_abs < th;
    r_causal = r(I_causal);
    p_causal = poly(r_causal);
    y = filtfilt(sum(p_causal), p_causal, x')';
elseif(mode == 1) % gamma-based
    gg = phi;
    gg(midindex) = gg(midindex) + SmoothnessFactor;
    r = roots(gg);
    r_abs = abs(r);
    I_causal = r_abs < th;
    r_causal = r(I_causal);
    p_causal = poly(r_causal);
    y = filtfilt(sum(p_causal), p_causal, x')'; % Note: the sqrt(SmoothnessFactor) which should be theoretically multiplied in x in this line is automatically compensated by the sum(p_causal) term
end
