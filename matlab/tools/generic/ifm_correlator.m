function f = ifm_correlator(x, f0, BW, tau)
% ifm_correlator - An instantaneous frequency estimator based on lagged correlations
%
% Syntax: f = ifm_correlator(x, f0, BW, tau)
%
% Inputs:
%   x: input signal as a row vector
%   f0: initial estimate of the frequency
%   BW: bandwidth for low-pass filtering
%   tau: time lag
%
% Output:
%   f: estimated instantaneous frequency
%
% References:
%   "Tsui, James, and Chi-Hao Cheng. Digital Techniques for Wideband
%   Receivers. 3rd edition, Scitech Publishing, 2016."
% 
%   Revision History:
%       2019: First release
%       2023: Renamed from deprecated version IFMCorrelator()
% 
%   Reza Sameni, 2019-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET


x = x(:)'; % Ensure x is a row vector

T = length(x);
n = (0:T-1);

c = cos(2*pi*f0*n); % Generate a cosine waveform
s = -sin(2*pi*f0*n); % Generate a negated sine waveform

y_i = x .* c; % Multiply input signal with the cosine waveform
y_q = x .* s; % Multiply input signal with the negated sine waveform

z_i = lp_filter_zero_phase(y_i, BW); % Apply low-pass filtering to the cosine component
z_q = lp_filter_zero_phase(y_q, BW); % Apply low-pass filtering to the sine component

z_i_lagged = [z_i(tau+1:end) zeros(1, tau)]; % Shift the filtered cosine component by tau samples
z_q_lagged = [z_q(tau+1:end) zeros(1, tau)]; % Shift the filtered sine component by tau samples

A1 = z_i_lagged - z_q; % Perform subtraction of the lagged cosine component and the sine component
B1 = z_q_lagged - z_i; % Perform subtraction of the lagged sine component and the cosine component
C1 = z_q_lagged + z_q; % Perform addition of the lagged sine component and the current sine component
D1 = -z_i_lagged + z_i; % Perform subtraction of the negated lagged cosine component and the current cosine component

A2 = A1.^2; % Square the result of A1
B2 = B1.^2; % Square the result of B1
C2 = C1.^2; % Square the result of C1
D2 = D1.^2; % Square the result of D1

A3 = lp_filter_zero_phase(A2, BW); % Apply low-pass filtering to A2
B3 = lp_filter_zero_phase(B2, BW); % Apply low-pass filtering to B2
C3 = lp_filter_zero_phase(C2, BW); % Apply low-pass filtering to C2
D3 = lp_filter_zero_phase(D2, BW); % Apply low-pass filtering to D2

E = B3 - A3; % Perform subtraction of B3 and A3
F = C3 - D3; % Perform subtraction of C3 and D3

f = atan2(E, F) / (2*pi*tau); % Calculate the estimated instantaneous frequency using the arctangent
