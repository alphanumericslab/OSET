function M = warping_transform(in_vector_boundaries, out_vector_boundaries, order)
% warping_transform - Generates matrices for time-warping one time-series to another (with different lengths).
% This function can be used for ECG phase wrapping, mapping ECG beats of different lengths
% into ECG beats of a fixed length, also known as the "phase domain" (see reference).
%
% Usage:
%   M = warping_transform(in_vector_boundaries, out_vector_boundaries, order)
%   Left multiplication of this matrix in a vector expands/compresses the
%   vector to the target length. The resulting vector is guaranteed to
%   match certain knot points between the input and output vectors.
%
% Inputs:
%   in_vector_boundaries: Scalar or monotonically increasing vector of
%       input vector knot point indices. When scalar, in_vector_boundaries
%       indicates the input length (the column space)
%   out_vector_boundaries: Scalar or or monotonically increasing vector of
%       output vector knot point indices. When scalar, out_vector_boundaries
%       indicates the output length (the row space)
%   order: Interpolation order (1 for linear, 2 for quadratic, 3 for 3rd order).
%
% Outputs:
%   M: Warping matrix for time-warping the input vector to match the output
%      vector's specified knot points.
%
% Note: The resulting vector (post multiplication by M) is guaranteed to
%   match certain knot points between the input and output vectors as a
%   result of multiplying this matrix. This includes the first, last, and
%   user-defined intermediate points listed in in_vector_boundaries and
%   out_vector_boundaries.
% 
% Example 1:
%   % Generate warping matrix using quadratic interpolation
%   in_len = 50;
%   out_len = 57;
%   order = 2;
%   M = warping_transform(in_len, out_len, order);
% 
% Example 2:
%   % Define input and output knot points
%   in_knots = [1, 5, 10, 15];
%   out_knots = [1, 7, 14, 20];
%
%   % Generate warping matrix using linear interpolation
%   order = 1;
%   M = warping_transform(in_knots, out_knots, order);
%   % The intermediate points are guaranteed to be mapped to one another
%
%   % Create an input vector
%   input_vector = cumsum(rand(1000, 1));
%
%   % Warp the input vector using the generated matrix
%   warped_vector = M * input_vector;
%
% 
% Method: The function generates a matrix of fractional sampling
%   coefficients (polynomial of order 'order').
% 
% Reference (usage):
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%   electrocardiogram decomposition using periodic component analysis. IEEE
%   Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2023: First release
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check if the number of input and output knot points match
if length(in_vector_boundaries) ~= length(out_vector_boundaries)
    error('input and output vectors should have the same number of knots');
end

% Check if the knot points start from 1
if in_vector_boundaries(1) ~= 1
    in_vector_boundaries = [1, in_vector_boundaries];
end

if out_vector_boundaries(1) ~= 1
    out_vector_boundaries = [1, out_vector_boundaries];
end

% Check for monotonicity of knot sequences
if ~isempty(find(diff(in_vector_boundaries) <= 0, 1))
    error('input knot sequence should be monotonically increasing')
end

if ~isempty(find(diff(out_vector_boundaries) <= 0, 1))
    error('output knot sequence should be monotonically increasing')
end

% Initialize the warping matrix
M = zeros(out_vector_boundaries(end), in_vector_boundaries(end));

% Iterate through knot points
for k = 1 : length(in_vector_boundaries) - 1
    in_len = in_vector_boundaries(k+1) - in_vector_boundaries(k) + 1;
    out_len = out_vector_boundaries(k+1) - out_vector_boundaries(k) + 1;
    fractional_update = (in_len -1) / (out_len - 1);

    % Select interpolation method based on 'order'
    switch order
        case 1 % linear interpolation
            t = in_vector_boundaries(k);
            for row = out_vector_boundaries(k) : out_vector_boundaries(k+1)
                col1 = floor(t);
                col2 = col1 + 1;
                if col1 <= in_vector_boundaries(end)
                    M(row, col1) = col2 - t ;
                end

                if col2 <= in_vector_boundaries(end)
                    M(row, col2) = t - col1;
                end
                t = t + fractional_update;
            end
        case 2 % quadratic interpolation
            t = in_vector_boundaries(k);
            for row = out_vector_boundaries(k) : out_vector_boundaries(k+1)
                col1 = ceil(t - 1.5 + eps);
                col2 = col1 + 1;
                col3 = col1 + 2;
                A = [col1.^2, col1, 1; col2.^2, col2, 1; col3.^2, col3, 1];
                coefs = [t.^2, t, 1] * pinv(A);
                if col1 >= 1 && col1 <= in_vector_boundaries(end)
                    M(row, col1) = coefs(1);
                end
                if col2 >= 1 && col2 <= in_vector_boundaries(end)
                    M(row, col2) = coefs(2);
                end
                if col3 >= 1 && col3 <= in_vector_boundaries(end)
                    M(row, col3) = coefs(3);
                end
                t = t + fractional_update;
            end
        case 3 % 3rd order interpolation
            t = in_vector_boundaries(k);
            for row = out_vector_boundaries(k) : out_vector_boundaries(k+1)
                col1 = ceil(t - 2.0 + eps);
                col2 = col1 + 1;
                col3 = col1 + 2;
                col4 = col1 + 3;
                A = [col1.^3, col1.^2, col1, 1; col2.^3, col2.^2, col2, 1; col3.^3, col3.^2, col3, 1; col4.^3, col4.^2, col4, 1];
                coefs = [t.^3, t.^2, t, 1] * pinv(A);
                if col1 >= 1 && col1 <= in_vector_boundaries(end)
                    M(row, col1) = coefs(1);
                end
                if col2 >= 1 && col2 <= in_vector_boundaries(end)
                    M(row, col2) = coefs(2);
                end
                if col3 >= 1 && col3 <= in_vector_boundaries(end)
                    M(row, col3) = coefs(3);
                end
                if col4 >= 1 && col4 <= in_vector_boundaries(end)
                    M(row, col4) = coefs(4);
                end
                t = t + fractional_update; 
            end
        otherwise
            error('Undefined warping model order');
    end
end