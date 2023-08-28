function [P, Q, R, S, T] = ecgWavesSoI_cRR(ecg, Rpeak, RRf)
%
% Approximates peaks and borders of the P, Q, R, S, and T waves in the input ECG.
%
% Usage:
%   [P, Q, R, S, T] = ecgWavesSoI_cRR(ecg, Rpeak, RRf)
%
% Inputs:
%   ecg: Single-channel ECG signal
%   Rpeak: R-peak positions (if RRf is false) or R-peak-referenced positions (if RRf is true)
%   RRf: R-peak referenced flag. If true, R peak is the reference for all the output points.
%      If false, they are based on the sample indices.
%
% Outputs:
%   P, Q, R, S, T: Each output is a matrix where the first, second, and third columns respectively
%       contain the onset, peak, and offset of the P, Q, R, S, T waves.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%       Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%       on Signal Processing 70 (2022): 3181-3192.
% 
% Revision History:
%   2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Ensure Rpeak is a column vector
if length(Rpeak) == length(ecg)
    Rpeak = find(Rpeak);
end

ecg = ecg(:); 
Rpeak = Rpeak(:);
RR = Rpeak(2:end) - Rpeak(1:end-1); 
RR = [median(RR); RR; median(RR)];

% Initialize arrays for each wave segment
P = nan(length(Rpeak), 3);
Q = nan(length(Rpeak), 3);
R = nan(length(Rpeak), 3);
S = nan(length(Rpeak), 3);
T = nan(length(Rpeak), 3);

% Iterate through each RR interval
for i = 1 : length(RR) - 1
    % P wave
    P(i, 1) = Rpeak(i) - floor(0.4 * RR(i));
    P(i, 3) = Rpeak(i) - floor(0.1 * RR(i));
    if P(i, 1) > 0
        [~, I] = max(abs(ecg(P(i, 1) : P(i, 3))));
        P(i, 2) = I + P(i, 1) - 1;
    end
    
    % Q wave
    Q(i, 1) = Rpeak(i) - floor(0.1 * RR(i));
    [~, I] = min(ecg(Rpeak(i) - floor(0.1 * RR(i)) : Rpeak(i)));
    Q(i, 2) = I + Q(i, 1) - 1;
    [~, I] = min(abs(ecg(Q(i, 2) : Rpeak(i))));
    Q(i, 3) = I + Q(i, 2);
    [~, I] = min(abs(ecg(Rpeak(i) - floor(0.1 * RR(i)) : Q(i, 2))));
    Q(i, 1) = I + Rpeak(i) - floor(0.1 * RR(i)) - 1;

    % S wave
    S(i, 3) = Rpeak(i) + floor(0.15 * (RR(i + 1)));
    S(i, 1) = Rpeak(i) + floor(0.02 * (RR(i + 1)));
    [~, I] = min(ecg(S(i, 1) : S(i, 3)));
    S(i, 2) = I + S(i, 1) - 1;

    % R wave
    R(i, 1) = Q(i, 3);
    R(i, 3) = S(i, 1);
    R(i, 2) = Rpeak(i);

    % T wave
    T(i, 3) = Rpeak(i) + floor(0.6 * RR(i + 1));
    T(i, 1) = Rpeak(i) + floor(0.15 * RR(i + 1));
    [~, I] = max(ecg(T(i, 1) : T(i, 3)));
    T(i, 2) = I + T(i, 1) - 1;
end

% Handle cases where wave segments go beyond signal bounds
P(P < 1 | P > length(ecg)) = nan;
Q(Q < 1 | Q > length(ecg)) = nan;
R(R < 1 | R > length(ecg)) = nan;
S(S < 1 | S > length(ecg)) = nan;
T(T < 1 | T > length(ecg)) = nan;

% Apply R-peak referencing if required
if RRf
    P = P - Rpeak;
    Q = Q - Rpeak;
    R = R - Rpeak;
    S = S - Rpeak;
    T = T - Rpeak;
end

end
