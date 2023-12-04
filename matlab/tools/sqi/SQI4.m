function [valindSQI, indSQI] = SQI4(x, Fs)
% 
% [valindSQI, indSQI] = SQI4(x, Fs)
% Channel selection based on spectral energy ratio
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   x: input JADE component of fetal ECG record
%   Fs: sampling Frequency (Hz)
% output:
%   valindSQI: values of "spectral energy ratio" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "spectral energy ratio" measure
%
% The Open Source Electrophysiological Toolbox, version 3.14, February 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/




% Comment out for individual using
% % %     Remove the mean
% % mn = mean(x,2)*ones(1,size(x,2));
% % x = x - mn;
% % 
% % %     Normalize the variance of x
% % x = x./var(x,0,2);

% Parameter nitializing
[CH, ~] = size(x);
% Fs = 1000;                                                  % Sampling Frequency (Hz)
Fn = Fs/2;                                                  % Nyquist Frequency (Hz)
Wp = [1.0   70]/Fn;                                         % Passband Frequency (Normalised)
Ws = [0.5   71]/Fn;                                         % Stopband Frequency (Normalised)
Rp =   1;                                                   % Passband Ripple (dB)
Rs = 150;                                                   % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                             % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                                  % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                                % Convert To Second-Order-Section For Stability
ind = zeros(1, CH);
yk = zeros(size(x));

for ch = 1: CH
% yk(ch, :) = filtfilt(sosbp, gbp*.7, x(ch, :)); %times .7 for gain
yk(ch, :) = filtfilt(sosbp, gbp, x(ch, :));                               % Filter Signal
ind(ch) = mean(yk(ch, :).^2)/mean(x(ch, :).^2);
end

% ranking channels according to the "spectral energy ratio" measure
[valindSQI, indSQI] = sort(ind, 'descend');