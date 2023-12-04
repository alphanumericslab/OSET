function [P, Q, R, S, T] = ecgWavesSoI_k(ecg, Rpeak, RRf)
%
% 'knot based' version.
% This function approximately finds the peaks and borders of the Q and T waves 
% in the input ECG signal.
%
% Inputs:
%   ecg: Single channel ECG signal
%   Rpeak: R peak positions
%   RRf: R-peak referenced flag; if true, R peak is the reference of all the
%       output points. If false, they are based on the sample index.
%
% Outputs:
%   P, Q, R, S, T: Each of the outputs is a matrix, where the first, second, and third 
%       columns respectively contain the onset, peak, and offset of the P, Q, R, S, and T waves.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Ensure Rpeak is in column vector form
if length(Rpeak) == length(ecg)
    Rpeak = find(Rpeak);
end

% Ensure ecg and Rpeak are column vectors
ecg = ecg(:);
Rpeak = Rpeak(:);

% Calculate RR intervals
RR = Rpeak(2:end) - Rpeak(1:end-1);
RR = [median(RR); RR; median(RR)];

% Initialize output matrices
P = nan(length(Rpeak), 3);
Q = nan(length(Rpeak), 3);
R = nan(length(Rpeak), 3);
S = nan(length(Rpeak), 3);
T = nan(length(Rpeak), 3);

% Loop through each RR interval
for i = 1:length(RR)-1
    % Q wave
    s = Rpeak(i) - floor(0.4 * RR(i));
    e = Rpeak(i);
    s = max(1, s);
    e = min(length(ecg), e);
    k = findknots(ecg(s:e), 3, 'Value') + s;
    Q(i, 2) = k(2);
    Q(i, 1) = Q(i, 2) - ceil(0.050 * RR(i));
    [~, I] = min(abs(ecg(Q(i, 2):Rpeak(i)) - ecg(Q(i, 1)))); % Nearest point to iso-point
    Q(i, 3) = Q(i, 2) + I - 1;
    
    % P wave
    s = Rpeak(i) - floor(0.4 * RR(i));
    e = Q(i, 1);
    s = max(1, s);
    e = min(length(ecg), e);
    k = findknots(ecg(s:e), 3, 'Value') + s;
    P(i, :) = k;
    
    % S wave
    s = Rpeak(i);
    e = Rpeak(i) + ceil(0.15 * RR(i + 1));
    s = max(1, s);
    e = min(length(ecg), e);
    k = findknots(ecg(s:e), 3, 'Value') + s;
    S(i, 2) = k(2);
    S(i, 3) = S(i, 2) + ceil(0.050 * RR(i));
    [~, I] = min(abs(ecg(s:S(i, 2)) - ecg(S(i, 3)))); % Nearest point to ST level
    S(i, 1) = s + I - 1;
    
    % R wave
    R(i, 1) = Q(i, 3);
    R(i, 2) = Rpeak(i);
    R(i, 3) = S(i, 1);
    
    % T wave
    s = Rpeak(i) + floor(0.1 * RR(i + 1));
    e = Rpeak(i) + floor(0.6 * RR(i + 1));
    s = max(1, s);
    e = min(length(ecg), e);
    k = findknots(ecg(s:e), 3, 'Value') + s;
    T(i, 2) = k(2);
    T(i, 3) = T(i, 2) + floor(0.25 * RR(i));
    T(i, 1) = T(i, 2) - floor(0.20 * RR(i));
    
    % 
end

% Set invalid points to NaN
P(P < 1 | P > length(ecg)) = nan;
Q(Q < 1 | Q > length(ecg)) = nan;
R(R < 1 | R > length(ecg)) = nan;
S(S < 1 | S > length(ecg)) = nan;
T(T < 1 | T > length(ecg)) = nan;

% If RRf is true, adjust outputs based on Rpeak
if RRf
    P = P - Rpeak;
    Q = Q - Rpeak;
    R = R - Rpeak;
    S = S - Rpeak;
    T = T - Rpeak;
end