function [x, fs, status] = ReadBinaryFile(ifname)
% ReadBinaryFile has been deprecated. Use read_binary_signal instead.
warning('ReadBinaryFile has been deprecated. Use read_binary_signal instead.');
[x, fs, status] = read_binary_signal(ifname);