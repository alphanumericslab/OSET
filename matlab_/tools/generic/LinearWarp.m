function y = LinearWarp(x, L)
% y = LinearWarp(x, L)
% Linear time warping of vectors/matrices to arbitrary lengths
%
% Reza Sameni
% Oct 2022
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET
%
 
if isvector(x)
    M = length(x);
    tx = (0 : M - 1) / (M - 1);
    ty = (0 : L - 1) / (L - 1);
    y = interp1(tx, x, ty);
elseif ismatrix(x)
    M1 = size(x, 1);
    M2 = size(x, 2);
    [t1, t2] = meshgrid((0 : M2 - 1) / (M2 - 1), (0 : M1 - 1) / (M1 - 1));
    [ty1, ty2] = meshgrid((0 : L(2) - 1) / (L(2) - 1), (0 : L(1) - 1) / (L(1) - 1));
    y = interp2(t1, t2, x, ty1, ty2);
else
    error('First input should be either a vector or a matrix');
end