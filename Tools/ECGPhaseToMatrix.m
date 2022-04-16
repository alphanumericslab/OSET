function [I, end_of_beat_indexes] = ECGPhaseToMatrix(phase, N)
% Converts ECG phase signal of length T into a matrix of N phase bins with the same length
% (used to calculate a binary image/matrix representation of the cardiac phase)
%
% Copyright Reza Sameni, 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

T = length(phase);
I = zeros(N, T);

phase_min = min(phase);
phase_max = max(phase);

loc = round((N - 1) * (phase - phase_min) / (phase_max - phase_min)) + 1;
for k = 1 : T
    I(loc(k), k) = 1;
end

end_of_beat_indexes = find(abs(diff(phase)) > pi);
% if end_of_beat_indexes(end) < T
%     end_of_beat_indexes = cat(2, end_of_beat_indexes, T);
% end
