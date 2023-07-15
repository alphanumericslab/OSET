function index1 = ChannelIndex13(x,ff,fs, wlen)
% Channel selection based on PiCA index
%
% Reza Sameni
% February 2019

L1 = size(x,1);
L2 = size(x,2);

index1 = zeros(L1,1);
mn = mean(x,2)*ones(1,L2);
x = x - mn;

for i = 1:L1
    %     yy = filtfilt(ones(1, wlen), wlen, x(i,:).^2);
    yy = x(i,:);
    peaks = peak_detection_local_search(yy, ff/fs, 1);
    [tt0, tt1] = SynchPhaseTimes2(peaks);

    % periodic component analysis stage
    A = x(:,tt0)*x(:,tt1)'/length(tt0);
    B = x(:,tt0)*x(:,tt0)'/length(tt0);

    A = (A+A')/2;
    B = (B+B')/2;

    index1(i) = trace(A)/trace(B);
end