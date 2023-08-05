function index = ChannelIndex13(x,ff,fs, wlen)
% ChannelIndex13 has been deprecated. Use sqi_pseudo_periodicity instead.
warning('ChannelIndex13 has been deprecated. Use sqi_pseudo_periodicity instead.');
method = 'TRACE';
num_peak_detection_itr = 1;
[index, ~] = sqi_pseudo_periodicity(x, ff, fs, method, num_peak_detection_itr);