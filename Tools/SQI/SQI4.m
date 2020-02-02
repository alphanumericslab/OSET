function ind = SQI4(x)
% Channel selection based on spectral energy ratio
%
% Fahimeh Jamshidian Tehrani
% February 2020

% Comment out for individual using
% % %     Remove the mean
% % mn = mean(x,2)*ones(1,size(x,2));
% % x = x - mn;
% % 
% % %     Normalize the variance of x
% % x = x./var(x,0,2);

[CH, ~] = size(x);


Fs = 1000;                                                  % Sampling Frequency (Hz)
Fn = Fs/2;                                                  % Nyquist Frequency (Hz)
Wp = [1.0   70]/Fn;                                         % Passband Frequency (Normalised)
Ws = [0.5   71]/Fn;                                         % Stopband Frequency (Normalised)
Rp =   1;                                                   % Passband Ripple (dB)
Rs = 150;                                                   % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                             % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                                  % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                                % Convert To Second-Order-Section For Stability

for ch = 1: CH
% yk(ch, :) = filtfilt(sosbp, gbp*.7, x(ch, :)); %times .7 for gain
yk(ch, :) = filtfilt(sosbp, gbp, x(ch, :));                               % Filter Signal
ind(ch) = mean(yk(ch, :).^2)/mean(x(ch, :).^2);
end
