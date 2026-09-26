
function [positions, EXITFLAG] = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)

% FIDUCIAL_DET_LSIM  (DEPRECATED) ECG fiducial points detector based on LSIM
%
%   [positions, EXITFLAG] = fiducial_det_lsim(data, ecg_rpeaks_index, fs, varargin)
%
%   This function is deprecated and kept only for backward compatibility.
%   It forwards all inputs to ECG_DELINEATE_LSIM and returns its
%   outputs. The optional inputs are identical, but ECG_DELINEATE_LSIM
%   takes fs as the second input:
%       fiducial_det_lsim(data, ecg_rpeaks_index, fs, ...)    % old
%       ecg_delineate_lsim(data, fs, ecg_rpeaks_index, ...)   % new
%   The output struct is a superset of the previous one (it additionally
%   contains beat_snr, index_clustering and rpeak_bp_lower_cutoff_hz).
%
%   Use ecg_delineate_lsim directly in new code. The previous
%   implementation is archived in depricated/fiducial_det_lsim_v2.m.
%
%   See also ECG_DELINEATE_LSIM.
%
%   Sajjad Karimi, Reza Sameni  2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
%

persistent warned
if isempty(warned)
    warning('OSET:fiducial_det_lsim:deprecated', ...
        ['fiducial_det_lsim is deprecated and now calls ecg_delineate_lsim. ' ...
        'Please call ecg_delineate_lsim directly.']);
    warned = true;
end

% note the input order of ecg_delineate_lsim: (data, fs, ecg_rpeaks_index, ...)
[positions, EXITFLAG] = ecg_delineate_lsim(data, fs, ecg_rpeaks_index, varargin{:});

end
