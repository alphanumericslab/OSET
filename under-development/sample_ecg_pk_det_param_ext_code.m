% Adding data, a single-channel ECG recording
f_hr = 1.2; % expected HR in Hz
fs = 1000; % Sampling frequency
num_rounds = 3; % number of times to iteratively run peak_det_local_search and tweak the peak locations

[~, peak_indexes] = peak_det_local_search(data, f_hr/fs, [], num_rounds); % peak_det_local_search is from OSET
heasig.nsig = 1; % number of channels
heasig.freq = fs; % 
heasig.nsamp = length(data);
[position, ~, ~] = wavedet_3D(data(:), peak_indexes, heasig); % wavedet_3D is from ecg-kit
