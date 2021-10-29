function I = ECGPhaseToMatrix(phase, N)
% Converts ECG phase signal of length T into a matrix of N phase bins with the same length
% (used to calculate a binary image/matrix representation of the cardiac phase)

% Reza Sameni, 2021
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
