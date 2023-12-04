function [A, AA, B] = PiCAMatrixes(x, peaks, peaks2)
% PiCAMatrixes has been deprecated. Use pica_matrices_double_peaks instead.
warning('PiCAMatrixes has been deprecated. Use pica_matrices_double_peaks instead.');
[A, AA, B] = pica_matrices_double_peaks(x, peaks, peaks2);