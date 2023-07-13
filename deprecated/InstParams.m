function [A, f0, bw, fskew, skew, kurt, f1] = InstParams(signal, fs, wlen, window, th, varargin)
% InstParams has been deprecated. Use instantaneous_signal_params instead.
warning('InstParams has been deprecated. Use instantaneous_signal_params instead.');
[A, f0, bw, fskew, skew, kurt, f1] = instantaneous_signal_params(signal, fs, wlen, window, th, varargin);