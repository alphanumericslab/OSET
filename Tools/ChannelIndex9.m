function ind = ChannelIndex9(x,ff,fs)
% Channel selection based on variance around average beat
%
% Reza Sameni
% November 2009

L1 = size(x,1);
L2 = size(x,2);
bins = round(0.25*fs);

ind = zeros(L1,1);
mn = mean(x,2)*ones(1,L2);
x = x - mn;

for i = 1:L1
    % % %     peaks = PeakDetection5(x(i,:),ff/fs);
    % % %     peaks = PeakDetectionFromC2(x(i,:), ff, fs);
    %     peaks = PeakDetectionFromC3(x(i,:), ff, fs);
    peaks = PeakDetection(x(i,:), ff/fs);
    % % %     phase = PhaseCalculation(peaks);
    % % %     [ECGmean,ECGsd] = MeanECGExtraction(x(i,:),phase,bins,1); % mean ECG extraction
    [ECGmean, ECGsd] = ECGBeatVariance(x(i,:), peaks, bins);
    ind(i) = max(abs(ECGmean))/mean(ECGsd);
    %     ind(i) = 1/mean(ECGsd);
    
    % % %     figure;
    % % %     plot(x(i,:),'k');
    % % %     grid;
    % % %     hold on;
    % % %     plot(peaks.*x(i,:),'ro');
    % % %     hold off;
end