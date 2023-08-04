function [x_filtered, x_filtered1] = ecg_den_seg_wise_smoother(x, varargin)
% ecg_den_seg_wise_smoother - Implementation of the piece-wise ECG filter using smoothness priors.
%
%   [x_filtered, x_filtered1] = ecg_den_seg_wise_smoother(x, DiffOrder, wlen, lambda, guardlen_l, guardlen_u, adapt, withwindow)
%
% Inputs:
%   x: Matrix of input signal (channels x samples).
%   DiffOrder: Integer representing the order of finite difference for the smoother. Default is 2.
%   wlen: Length of each segment used in the smoother. Default is approximately 10% of the signal length.
%   lambda: Regularization parameter controlling the smoothness of the smoother. Default is 1.0.
%   guardlen_l: Number of guard samples to the left of each segment during the second round. Default is 3.
%   guardlen_u: Number of guard samples to the right of each segment during the second round. Default is 3.
%   adapt: Boolean flag indicating whether to use an adaptive lambda based on segment length and data. Default is false.
%   withwindow: Boolean flag indicating whether to apply a Hamming window on each segment. Default is false.
%
% Outputs:
%   x_filtered: Signal matrix after the second round of filtering.
%   x_filtered1: Signal matrix after the first round of filtering (has discontinuities across segment borders).
%
% Note: This implementation corresponds to the fixed lambda implementation
%   of the piece-wise signal smoother from the reference below.
%
% Reference: See Section 2 of
%   R. Sameni, Online filtering using piecewise smoothness priors: Application
%   to normal and abnormal electrocardiogram denoising, Signal Processing,
%   Volume 133, 2017, Pages 52-63, https://doi.org/10.1016/j.sigpro.2016.10.019.
%
% Revision History:
%   2015: First release
%   2023: Renamed from deprecated version ECGSmoothingPriors and added
%       multichannel feature
%
% Reza Sameni, 2015-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Handle optional input arguments
if nargin > 1 && ~isempty(varargin{1})
    DiffOrder = varargin{1};
else
    DiffOrder = 2;
    warning(['Unspecified smoothness order set to: ', num2str(DiffOrder)]);
end

if nargin > 2 && ~isempty(varargin{2})
    wlen = varargin{2};
else
    wlen = ceil(size(x, 2)/10.0);
    warning(['Unspecified wlen set to: ', num2str(wlen)]);
end

if nargin > 3 && ~isempty(varargin{3})
    lambda = varargin{3};
else
    lambda = 1.0;
    warning(['Unspecified lambda set to: ', lambda]);
end

if nargin > 4 && ~isempty(varargin{4})
    guardlen_l = varargin{4};
else
    guardlen_l = 3;
end

if nargin > 5 && ~isempty(varargin{5})
    guardlen_u = varargin{5};
else
    guardlen_u = 3;
end

if nargin > 6 && ~isempty(varargin{6})
    adapt = varargin{6};
else
    adapt = false;
end

if nargin > 7 && ~isempty(varargin{7})
    withwindow = varargin{7};
else
    withwindow = false;
end

N = size(x, 1); % Number of channels
T = size(x, 2); % Number of samples

x_filtered1 = zeros(N, T); % Initialize the matrix to store the results after the first round
x_filtered = zeros(N, T); % Initialize the matrix to store the results after the second round

knots1 = 1 : round(wlen) : T; % Generate the starting indices of each segment

% Add the first and last samples as knots to avoid end-point discontinuities
if knots1(1) > 1
    knots1 = cat(2, 1, knots1);
end
if knots1(end) < T
    knots1 = cat(2, knots1, T);
end

knots2 = round(wlen / 2) : round(wlen) : T; % Generate the middle indices of each segment

if knots2(1) > 1
    knots2 = cat(2, 1, knots2);
end
if knots2(end) < T
    knots2 = cat(2, knots2, T);
end

for ch = 1:N % Loop through each channel

    xx = x(ch, :)'; % Extract the data for the current channel

    % First round of filtering
    for n = 1:length(knots1) - 1
        indexes = knots1(n):knots1(n + 1); % Get the current segment
        indexes_len = length(indexes);
        I = speye(indexes_len); % Identity matrix for the current segment
        % Create the finite difference matrix (Dd) to apply the smoother
        Dd = spdiags(ones(indexes_len - 2, 1) * diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, indexes_len - DiffOrder, indexes_len);
        %     R = qr(I + lambda^2*(Dd'*Dd));
        %     x_filtered1(indexes) = R\(R'\xx(indexes));

        % Apply the smoother using regularization parameter lambda
        if withwindow
            W = diag(hamming(length(I))); % Apply a Hamming window on each segment
        else
            W = I;
        end
        x_filtered1(ch, indexes) = (W + lambda^2 * (Dd' * Dd)) \ (W * xx(indexes));
    end

    % Second round of filtering
    for n = 1:length(knots1) - 1
        indexes = knots2(n):knots2(n + 1); % Get the current segment
        indexes_len = length(indexes);
        I = speye(indexes_len); % Identity matrix for the current segment size

        % Create the finite difference matrix (Dd) to apply the smoother
        Dd = spdiags(ones(indexes_len - 2 + guardlen_l + guardlen_u, 1) * diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, indexes_len + guardlen_l + guardlen_u - DiffOrder, indexes_len + guardlen_l + guardlen_u);
        D_l = Dd(:, 1:guardlen_l);
        D = Dd(:, guardlen_l + 1:end - guardlen_u);
        D_u = Dd(:, end - guardlen_u + 1:end);

        % Indices for the guard samples to the left and right of the current segment
        indxl = indexes(1) - guardlen_l : indexes(1) - 1;
        indxl(indxl < 1) = 1;

        indxu = indexes(end) + 1 : indexes(end) + guardlen_u;
        indxu(indxu > length(x_filtered1)) = length(x_filtered1);

        % Extract the guard samples
        l = x_filtered1(ch, indxl)';
        u = x_filtered1(ch, indxu)';

        % Calculate the adaptive lambda for the current segment
        if adapt
            lambda_adaptive = (lambda / (sqrt(T / indexes_len) * norm(xx(indexes)) / norm(xx))); % Adaptive version
        else
            lambda_adaptive = lambda; % Fixed version
        end

        % Apply the smoother using the adaptive lambda
        if withwindow
            W = diag(hamming(length(I))); % Apply a Hamming window on each segment
        else
            W = I;
        end
        %     R = qr(I + lambda_adaptive^2*(D'*D));
        %     x_filtered(indexes) = R\(R'\(xx(indexes) - lambda_adaptive^2*(D'*D_l(:, 1:length(l)))*l - lambda_adaptive^2*(D'*D_u(:, 1:length(u)))*u));
        x_filtered(ch, indexes) = (W + lambda_adaptive^2 * (D' * D)) \ ...
            (W * xx(indexes) - lambda_adaptive^2 * (D' * D_l(:, 1:length(l))) * l - lambda_adaptive^2 * (D' * D_u(:, 1:length(u))) * u);
    end

end