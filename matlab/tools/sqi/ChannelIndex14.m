function ind = ChannelIndex14(x,ff,fs, varargin)
% Channel selection based on variance around average beat
% same as ChannelIndex9(); but using ECGBeatVariance2() which makes the
% R-peaks equal in amplitude
%
% Reza Sameni
% November 2009

L1 = size(x,1);
L2 = size(x,2);
if(nargin == 3)
    bins = round(0.25*fs);
else
    bins = varargin{1};
end

ind = zeros(L1,1);
mn = mean(x,2)*ones(1,L2);
x = x - mn;

x = x./(std(x, [], 2)*ones(1, size(x, 2)));
for i = 1:L1
    % % %     peaks = PeakDetection5(x(i,:),ff/fs);
    % % %     peaks = PeakDetectionFromC2(x(i,:), ff, fs);
    %     peaks = PeakDetectionFromC3(x(i,:), ff, fs);
    peaks = PeakDetection(x(i,:), ff/fs);
    % % %     phase = PhaseCalculation(peaks);
    % % %     [ECGmean,ECGsd] = MeanECGExtraction(x(i,:),phase,bins,1); % mean ECG extraction
    [ECGmean, ECGsd] = ECGBeatVariance2(x(i,:), peaks, bins);
    ind(i) = max(abs(ECGmean))/mean(ECGsd);
    %     ind(i) = 1/mean(ECGsd);
    
    % % %     figure;
    % % %     plot(x(i,:),'k');
    % % %     grid;
    % % %     hold on;
    % % %     plot(peaks.*x(i,:),'ro');
    % % %     hold off;
end