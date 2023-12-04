function y = linear_warp(x, L)
% linear_warp - Linear time warping of vectors/matrices to arbitrary lengths
%   y = linear_warp(x, L)
%   Linear time warping of vectors/matrices to arbitrary lengths.
%
% Inputs:
%   x: input vector or matrix
%   L: length or size of the output vector/matrix
%
% Output:
%   y: time-warped output vector/matrix
%
% Note:
% - This function performs linear time warping of the input signal
%   to match the desired length or size specified by L.
% - If x is a vector, linear interpolation is used for time warping.
% - If x is a matrix, bilinear interpolation is used for time warping.
% - In the vector case, L is the desired length of the output vector.
% - In the matrix case, L is a 2-element vector specifying the desired
%   size of the output matrix as [num_rows, num_columns].
%
%   Revision History:
%       2022: First release
%       2023: Renamed from deprecated version LinearWarp()
% 
%   Reza Sameni, 2022-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
 
if isvector(x)
    % Time warping for vectors
    M = length(x);
    tx = (0 : M - 1) / (M - 1); % Time axis for the original signal
    ty = (0 : L - 1) / (L - 1); % Time axis for the desired output signal
    y = interp1(tx, x, ty); % Linear interpolation to warp the signal
elseif ismatrix(x)
    % Time warping for matrices
    M1 = size(x, 1); % Number of rows in the original matrix
    M2 = size(x, 2); % Number of columns in the original matrix
    
    % Interpretation of L in the matrix case
    num_rows = L(1); % Desired number of rows for the output matrix
    num_columns = L(2); % Desired number of columns for the output matrix
    
    [t1, t2] = meshgrid((0 : M2 - 1) / (M2 - 1), (0 : M1 - 1) / (M1 - 1)); % Time grid for the original matrix
    [ty1, ty2] = meshgrid((0 : num_columns - 1) / (num_columns - 1), (0 : num_rows - 1) / (num_rows - 1)); % Time grid for the desired output matrix
    y = interp2(t1, t2, x, ty1, ty2); % Bilinear interpolation to warp the matrix
else
    error('First input should be either a vector or a matrix');
end
