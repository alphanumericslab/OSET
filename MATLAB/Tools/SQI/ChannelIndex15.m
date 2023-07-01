function ind = ChannelIndex15(x, ff, fs, wlen, guard, th)
% Channel selection based on periodicity of local energies
%
% Reza Sameni
% February 2019

N = size(x,1);
T = size(x,2);
nfft = T;
ind = zeros(N, 1);
for l = 1 : N
    envelope = sqrt(filtfilt(ones(1, wlen), wlen, x(l,:).^2));
    
    %     envelope_pulse = zeros(1, T);
    %     indexes = envelope > th*max(envelope);
    %     envelope_pulse(indexes) = 1;
    
    ppkkss = PeakDetection(envelope, ff/fs, 1);
    %     ppkkss = PeakDetection6(envelope, ff/fs, th, 1);
    
    ff_updated = fs/median(diff(find(ppkkss)));
    
    XX = abs(fft(envelope, nfft));
    %     XX = abs(fft(envelope_pulse, nfft));
    kk0 = round(nfft*ff_updated/fs);
    kk = [kk0-guard:kk0+guard 2*kk0-guard:2*kk0+guard 3*kk0-guard:3*kk0+guard 4*kk0-guard:4*kk0+guard 5*kk0-guard:5*kk0+guard]; % find the first few harmonics of the peak
    kk(kk > floor(nfft/2)) = [];
    kk(kk < 1) = [];
    ind(l) = mean(XX(kk))/mean(XX(1:floor(nfft/2)));
end
