function [P, Q, R, S, T] = ecgWavesSoI_RRT(ecg, Rpeak, fs, RRf)
%
% 'RR interval & time based' version.
% This function approximately finds the peaks and borders of the Q and T waves 
% in the input ECG signal using RR intervals and time-based information.
%
% Inputs:
%   ecg: Single channel ECG signal
%   Rpeak: R peak positions
%   fs: Sampling frequency
%   RRf: R-peak referenced flag; if true, R peak is the reference of all the
%   output points. If false, they are based on the sample index.
%
% Outputs:
%   P, Q, R, S, T: Each of the outputs is a matrix, where the first, second, and third 
%   columns respectively contain the onset, peak, and offset of the P, Q, R, S, and T waves.
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
    % P wave
    P(i, 1) = Rpeak(i) - floor(0.4 * RR(i));
    P(i, 3) = Rpeak(i) - floor(0.1 * RR(i));
    [~, I] = max(abs(ecg(P(i, 1):P(i, 3))));
    P(i, 2) = I + P(i, 1) - 1;
    
    % Q wave
    Q(i, 1) = Rpeak(i) - floor(0.08 * fs);
    Q(i, 3) = Rpeak(i) - floor(0.02 * fs);
    [~, I] = min(ecg(Q(i, 1):Q(i, 3)));
    Q(i, 2) = I + Q(i, 1) - 1;
    
    % S wave
    S(i, 3) = Rpeak(i) + floor(0.15 * RR(i + 1));
    S(i, 1) = Rpeak(i) + floor(0.04 * RR(i + 1));
    [~, I] = min(ecg(S(i, 1):S(i, 3)));
    S(i, 2) = I + S(i, 1) - 1;
    
    % R wave
    R(i, 1) = Q(i, 3);
    R(i, 3) = S(i, 1);
    R(i, 2) = Rpeak(i);
    
    % T wave
    T(i, 3) = Rpeak(i) + floor(0.6 * RR(i + 1));
    T(i, 1) = Rpeak(i) + floor(0.1 * fs);
    [~, I] = max(ecg(T(i, 1):T(i, 3)));
    T(i, 2) = I + T(i, 1) - 1;
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

% Uncomment the following lines to plot the detected points
% fidp = [Q(:); R(:); S(:); T(:)];
% figure;
% plot(ecg); hold on; plot(fidp, ecg(fidp), 'r*');

end
