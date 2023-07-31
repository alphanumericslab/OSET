clear;
clc;
close all;

fs = 1200.0;
x = randn(5, 500);
fname = 'sample_binary_file.bin';
[status1, count] = MatToBinaryFileConverter(x, fs, fname);

[xx, ffs, status2] = ReadBinaryFile(fname);

isequal(fs, ffs)

isequal(x, xx)
