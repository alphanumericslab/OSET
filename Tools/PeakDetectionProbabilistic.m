function [peaks, I_peaks, qrs_likelihood] = PeakDetectionProbabilistic(signal, fs, varargin)
% A probabilistic R-peak detector based on local peaks sorting
%
% Usage:
%   [peaks, qrs_likelihood, I_peaks] = PeakDetectionProbabilistic(signal, fs, params)
%
% Inputs:
%   signal: single or multichannel ECG signal with row-wise channels
%   fs: sampling frequency
%   params: a structure containing the algorithm parameters (all parameters have default values if not given as input)
%       params.saturate: saturate (1) the channels before processing or not (0). Default = 0
%       params.low_cutoff: lower cutoff frequency of R-peak detector bandpass filter (in Hz). Default = 15
%       params.up_cutoff: upper cutoff frequency of R-peak detector bandpass filter (in Hz). Default = 40
%       params.wlen_power_env: power envelope moving average window length (in seconds). Default = 0.03
%       params.n_bins: number of local peaks histogram bins. Default = max(250, min(500, 10% of signal length))
%       params.hist_search_th: the top percentile of the signal's power envelope, considered for R-peak detection. Default = 0.9 (top 10%)
%       params.rpeak_search_wlen: the R-peak detector search window length (in seconds). Default = 0.1
%       params.likelihood_sigma: R-peak likelihood estimator STD in seconds (assuming a Gaussian that drops from the peak maximum)
%       params.max_likelihood_span: R-peak likelihood estimator max window span in seconds (assuming a Gaussian that drops from the peak maximum)
%
% Outputs:
%   peaks: a vector with the signal length with 1's at the estimated R-peaks
%   I_peaks: the estimated R-peak indexes
%   qrs_likelihood: the R-peak likelihood vector (with maximums at the
%       estimated R-peaks useful for classification and scoring purposes)
%
% Copyright Reza Sameni, Nov 2021
% The Open-Source Electrophysiological Toolbox
% (https://github.com/alphanumericslab/OSET)

% use default values when no parameters are set
if nargin > 2
    params = varargin{1};
else
    params = [];
end

SigLen = size(signal, 2); % signal length

% Saturate the channels at k_sigma times the channel-wise STD
if isfield(params, 'saturate') && isequal(params.saturate, 1)
    if ~isfield(params, 'k_sigma')
        params.k_sigma = 8.0;
    end
    alpha = k_sigma * std(signal, [], 2);
    alpha = alpha(:);
    data_sat = diag(alpha) * tanh(diag(1./alpha) * signal);
else
    data_sat = signal;
end


% pass the channels through a narrow bandpass filter
if ~isfield(params, 'low_cutoff')
    params.low_cutoff = 15.0;
end
if ~isfield(params, 'up_cutoff')
    params.up_cutoff = 40.0;
end
data_bp = LPFilter(data_sat - LPFilter(data_sat, params.low_cutoff/fs), params.up_cutoff/fs);

% calculate the power envelope of one or all channels
if ~isfield(params, 'wlen_power_env')
    params.wlen_power_env = 0.03;
end
wlen_power_env = ceil(params.wlen_power_env * fs);
data_bp_env = filtfilt(ones(1, wlen_power_env), wlen_power_env, sqrt(mean(data_bp.^2, 1)));

% calculate the signal's power envelope histogram
if ~isfield(params, 'n_bins')
    params.n_bins = min(500, max(250, round(SigLen/10))); % somewhere between 250 and 500, if not specified
end
[data_bp_env_hist.Values, data_bp_env_hist.BinEdges] = histcounts(data_bp_env, params.n_bins);
data_bp_env_PDF = data_bp_env_hist.Values/sum(data_bp_env_hist.Values); % calculate the PDF
data_bp_env_CDF = cumsum(data_bp_env_PDF); % calculate the CDF

% pick the top percentage of the signal's power envelope for R-peak search
if ~isfield(params, 'hist_search_th')
    params.hist_search_th = 0.9;
end
bumps_threshold_index = find(data_bp_env_CDF >= params.hist_search_th, 1, 'first');
bumps_amp_threshold = data_bp_env_hist.BinEdges(bumps_threshold_index);
bumps_indexes = find(data_bp_env >= bumps_amp_threshold);

% search among the top percentile for the local peaks within a given sliding window length
if ~isfield(params, 'rpeak_search_wlen')
    params.rpeak_search_wlen = 0.1;
end
rpeak_search_half_wlen = floor(fs * params.rpeak_search_wlen / 2);
I_peaks = [];
for jj = 1 : length(bumps_indexes)
    index = bumps_indexes(jj);
    segment = max(1, index - rpeak_search_half_wlen) : min(SigLen, index + rpeak_search_half_wlen);
    if max(data_bp_env(segment)) == data_bp_env(index)
        I_peaks = cat(2, I_peaks, index);
    end
end
peaks = zeros(1, SigLen);
peaks(I_peaks) = 1;

% calculate a likelihood function for the R-peaks (useful for classification and scoring purposes)
if ~isfield(params, 'likelihood_sigma')
    params.likelihood_sigma = 0.01;
end
if ~isfield(params, 'max_likelihood_span')
    params.max_likelihood_span = 0.4;
end
tt = -params.max_likelihood_span/2 : 1 / fs : params.max_likelihood_span/2;
template = exp(-tt.^2/(2 * params.likelihood_sigma ^ 2));
lag = round(length(template)/2);
qrs_likelihood = conv(template, peaks);
qrs_likelihood = qrs_likelihood(lag : SigLen + lag - 1);
