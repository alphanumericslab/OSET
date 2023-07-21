function [yy, criterion] = PeriodicDeflDecompositionWDEN(dat, Itr, MM, T0, T1, TPTR, SORH, SCAL, NDEN, WNAME)
%
% [yy criterion] = PeriodicDeflDecompositionWDEN(dat, Itr, MM, T0, T1, TPTR, SORH, SCAL, NDEN, WNAME)
% (maternal) cardiac signal supression by deflation using periodic component
% analysis and wavelet shrinkage for nonlinear denoising.
%
% Ref: An implementation of the following paper:
% Sameni, Reza, Christian Jutten, and Mohammad B. Shamsollahi.
% "A deflation procedure for subspace decomposition."
% IEEE Transactions on Signal Processing 58, no. 4 (2010): 2363-2374.
%
% Open Source ECG Toolbox, version 3.14, November 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

L2 = size(dat,2);

% mean removal
mn = mean(dat,2)*ones(1,L2);
dat = dat - mn;
criterion = zeros(1, Itr);
for i = 1 : Itr

    % periodic component analysis stage
    A = dat(:,T0)*dat(:,T1)'/length(T0);
    B = dat(:,T0)*dat(:,T0)'/length(T0);

    A = (A+A')/2;
    B = (B+B')/2;

    criterion(i) = trace(A)/trace(B);

    [V,D] = eig(A,B,'chol');

    d = diag(D);
    [~,I] = sort(d);
    I = I(end:-1:1);

    W = V(:,I)';
    dat = W*dat;

    % wavelet denoising
    for k = 1:MM
        est = wden(dat(k,:), TPTR, SORH, SCAL, NDEN, WNAME);
        dat(k,:) = dat(k,:) - est;
    end

    % reconstruction
    dat = pinv(W)*dat;
end

% return data mean
dat = dat + mn;

yy = dat;