function varargout = ECGBeatFitter(varargin)
%
% ECGBeatFitter(ECG,Phase,ExpParamName),
% Graphical user interface for ECG approximation with Gaussian kernels.
% 
% ECGBeatFitter has been deprecated. Use ecg_beat_fitter_gmm_gui instead.
warning('ECGBeatFitter has been deprecated. Use ecg_beat_fitter_gmm_gui instead.');

varargout = ecg_beat_fitter_gmm_gui(varargin{:});
