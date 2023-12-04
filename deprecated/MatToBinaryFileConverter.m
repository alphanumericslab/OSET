function [status, count] = MatToBinaryFileConverter(x, fs, ofname)
% MatToBinaryFileConverter has been deprecated. Use write_binary_signal instead.
warning('MatToBinaryFileConverter has been deprecated. Use write_binary_signal instead.');
[status, count] = write_binary_signal(x, fs, ofname);