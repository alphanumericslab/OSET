clear
close all
clc

load ECG2
s = data(4, 1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'MATLAB:nearlySingularMatrix'); % THIS IS JUST FOR AVOIDING THE MATRIX SINGULARITY WARNINGS. A MORE NUMERICALLY STABLE VERSION OF THE CODE MIGHT BE NEEDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smoothingwlen1 = 3;
smoothingwlen2 = 3;
smoothingwlen3 = 3;
wlen1 = 150e-3;
wlen2 = 150e-3;
wlen3 = 150e-3;
p1 = 7; % spline order
p2 = 7; % spline order
p3 = 7; % spline order
snr = 10; % dB
f0 = 1.0; % Heart Beart in Hz
lmbd = 10;
gm = 50;

% th = 1e-45; % SVD inversion threshold

s = LPFilter(s - LPFilter(s, 0.3/fs), 60.0/fs);

nvar = var(s)/10^(snr/10);
x = s + sqrt(nvar)*randn(size(s));

% first round
II = 1: round(wlen1*fs) : length(x);
if(II(1) > 1)
    II = [1 II];
end
if(II(end) < length(x))
    II = [II length(x)];
end
y = zeros(size(x));
N = length(II);
for n = 1 : N - 1,
    i = II(n) : II(n+1); % I used this line to preserve the continuity of the curves
    ti = (0:length(i)-1)/fs;
    %     ti = (0:length(i)-1); % THIS AVOIDED THE SINGULARITIES
    t0 = [1 ; zeros(p1, 1)];
    tn = zeros(p1+1, 1);
    t = zeros(p1+1, length(i));
    for k = 1 : p1+1,
        tn(k) = ti(end)^(k-1);
        t(k, :) = ti.^(k-1);
    end
    
    Phi = [t0 tn];
    T = t*t';
    %     disp(['size(T) = ', num2str(length(T)), ', rank(T) = ', num2str(rank(T))]);
    %     S = Phi'*inv(T)*Phi;
    %     [U D V] = svd(T);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     T_inv = V*diag(d_inv)*U';
    
    S = Phi'/T*Phi;
    %     S = Phi'*T_inv*Phi;
    u = t*x(i)';
    
    %     [U D V] = svd(S);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     S_inv = V*diag(d_inv)*U';
    
    %     gamma = [x(i(1)) x(i(end))]';
    % this is an intutive tweak, to force the curve pass an average point
    % rather than the exact signal points at the curve ends:
    alpha = mean(x(max(i(1)-smoothingwlen1,1):min(i(1)+smoothingwlen1, length(x))));
    beta = mean(x(max(i(end)-smoothingwlen1,1):min(i(end)+smoothingwlen1, length(x))));
    gamma = [alpha ; beta];
    
    %     lambda = inv(S)*(Phi'*inv(T)*u - gamma);
    %     aa = inv(T)*(u - Phi*lambda);
    lambda = S\(Phi'/T*u - gamma);
    aa = T\(u - Phi*lambda);
    %     lambda = S_inv*(Phi'*T_inv*u - gamma);
    %     aa = T_inv*(u - Phi*lambda);
    y(i) = aa'*t;
end
% y = LPFilter(y, 80/fs);

% second round
II = round(wlen2*fs/2): round(wlen2*fs) : length(y);
if(II(1) > 1)
    II = [1 II];
end
if(II(end) < length(y))
    II = [II length(y)];
end
z = zeros(size(y));
N = length(II);
x_detrended0 = zeros(size(z));
for n = 1 : N - 1,
    i = II(n) : II(n+1); % I used this line to preserve the continuity of the curves
    ti = (0:length(i)-1)/fs;
    %     ti = (0:length(i)-1); % THIS AVOIDED THE SINGULARITIES
    t0 = [1 ; zeros(p2, 1)];
    t0_prime = [0 ; 1 ; zeros(p2 - 1, 1)];
    tn = zeros(p2+1, 1);
    tn_prime = zeros(p2+1, 1);
    t = zeros(p2+1, length(i));
    for k = 1 : p2+1,
        tn(k) = ti(end)^(k-1);
        tn_prime(k) = (k-1)*ti(end)^(k-2);
        t(k, :) = ti.^(k-1);
    end
    
    Phi = [t0 tn t0_prime tn_prime];
    T = t*t';
    
    %     [U D V] = svd(T);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     T_inv = V*diag(d_inv)*U';
    
    %     S = Phi'*inv(T)*Phi;
    S = Phi'/T*Phi;
    %     S = Phi'*T_inv*Phi;
    u = t*x(i)'; % the second stage is still applied to the original data; but we get the boundary conditions from the first stage
    
    %     [U D V] = svd(S);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     S_inv = V*diag(d_inv)*U';
    
    %     gamma = [x(i(1)) x(i(end))]';
    % this is an intutive tweak, to force the curve pass an average point
    % rather than the exact signal points at the curve ends:
    alpha = mean(y(max(i(1)-smoothingwlen2,1) : min(i(1)+smoothingwlen2, length(y))));
    beta = mean(y(max(i(end)-smoothingwlen2,1) : min(i(end)+smoothingwlen2, length(y))));
    alpha_prime = mean(fs*diff(y(max(i(1)-smoothingwlen2,1) : min(i(1)+smoothingwlen2, length(y)))));
    beta_prime = mean(fs*diff(y(max(i(end)-smoothingwlen2,1) : min(i(end)+smoothingwlen2, length(y)))));
    %     alpha_prime = mean(diff(y(max(i(1)-smoothingwlen2,1) : min(i(1)+smoothingwlen2, length(y)))));
    %     beta_prime = mean(diff(y(max(i(end)-smoothingwlen2,1) : min(i(end)+smoothingwlen2, length(y)))));
    gamma = [alpha ; beta ; alpha_prime ; beta_prime];
    
    %     lambda = inv(S)*(Phi'*inv(T)*u - gamma);
    %     aa = inv(T)*(u - Phi*lambda);
    lambda = S\(Phi'/T*u - gamma);
    aa = T\(u - Phi*lambda);
    %     lambda = S_inv*(Phi'*T_inv*u - gamma);
    %     aa = T_inv*(u - Phi*lambda);
    z(i) = aa'*t;

    ll = length(i);
    I = speye(ll);
    D2 = spdiags(ones(ll-2, 1)*[1 -2 1], 0:2, ll-2, ll);
    x_detrended0(i) = inv(I + lmbd^2*(D2'*D2))*x(i)';
end
% z = LPFilter(z, 80/fs);

% third round
II = 1: round(wlen3*fs) : length(z);
if(II(1) > 1)
    II = [1 II];
end
if(II(end) < length(z))
    II = [II length(z)];
end
w = zeros(size(z));
x_detrended = zeros(size(z));
x_detrended_mod = zeros(size(z));
N = length(II);
for n = 1 : N - 1,
    i = II(n) : II(n+1); % I used this line to preserve the continuity of the curves
    ti = (0:length(i)-1)/fs;
    %     ti = (0:length(i)-1); % THIS AVOIDED THE SINGULARITIES
    t0 = [1 ; zeros(p3, 1)];
    t0_prime = [0 ; 1 ; zeros(p3 - 1, 1)];
    t0_double_prime = [0 ; 0 ; 2 ; zeros(p3 - 2, 1)];
    tn = zeros(p3+1, 1);
    tn_prime = zeros(p3+1, 1);
    tn_double_prime = zeros(p3+1, 1);
    t = zeros(p3+1, length(i));
    for k = 1 : p3+1,
        tn(k) = ti(end)^(k-1);
        tn_prime(k) = (k-1)*ti(end)^(k-2);
        tn_double_prime(k) = (k-1)*(k-2)*ti(end)^(k-3);
        t(k, :) = ti.^(k-1);
    end
    
    Phi = [t0 tn t0_prime tn_prime t0_double_prime tn_double_prime];
    T = t*t';
    %     S = Phi'*inv(T)*Phi;
    
    %     [U D V] = svd(T);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     T_inv = V*diag(d_inv)*U';
    
    S = Phi'/T*Phi;
    %     S = Phi'*T_inv*Phi;
    u = t*x(i)'; % the third stage is still applied to the original data; but we get the boundary conditions from the first stage
    
    %     [U D V] = svd(S);
    %     d = diag(D);
    %     ii = find(d >= th);
    %     d_inv = zeros(size(d));
    %     d_inv(ii) = 1./d(ii);
    %     S_inv = V*diag(d_inv)*U';
    
    %     gamma = [x(i(1)) x(i(end))]';
    % this is an intutive tweak, to force the curve pass an average point
    % rather than the exact signal points at the curve ends:
    alpha = mean(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))));
    beta = mean(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))));
    alpha_prime = mean(fs*diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z)))));
    beta_prime = mean(fs*diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z)))));
    alpha_double_prime = mean(fs^2*diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))), 2));
    beta_double_prime = mean(fs^2*diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))), 2));
    %     alpha_prime = mean(diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z)))));
    %     beta_prime = mean(diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z)))));
    %     alpha_double_prime = mean(diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))), 2));
    %     beta_double_prime = mean(diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))), 2));
    gamma = [alpha ; beta ; alpha_prime ; beta_prime ; alpha_double_prime ; beta_double_prime];
    
    %     lambda = inv(S)*(Phi'*inv(T)*u - gamma);
    %     aa = inv(T)*(u - Phi*lambda);
    lambda = S\(Phi'/T*u - gamma);
    aa = T\(u - Phi*lambda);
    %     lambda = S_inv*(Phi'*T_inv*u - gamma);
    %     aa = T_inv*(u - Phi*lambda);
    w(i) = aa'*t;
    
    
    ll = length(i);
    I = speye(ll);
    D2 = spdiags(ones(ll-2, 1)*[1 -2 1], 0:2, ll-2, ll);
    x_detrended(i) = inv(I + lmbd^2*(D2'*D2))*x(i)';
    
    % my modified detrending that gives penalty to the start/end discontinuties
    G = zeros(6, ll);
    G(1,1) = 1;
    G(2,end) = 1;
    G(3, 1:2) = [-1 1];
    G(4, end-1:end) = [-1 1];
    G(5, 1:3) = [1 -2 1];
    G(6, end-2:end) = [1 -2 1];
    alpha = x_detrended0(i(1));%mean((max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))));
    beta = x_detrended0(i(end));%mean(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))));
    alpha_prime = diff(x_detrended0(i(1:2)));%mean(diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z)))));
    beta_prime = diff(x_detrended0(i(end-1:end)));%mean(diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z)))));
    alpha_double_prime = diff(x_detrended0(i(1:3)), 2);%mean(diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))), 2));
    beta_double_prime = diff(x_detrended0(i(end-2:end)), 2);%mean(diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))), 2));
    gamma = [alpha ; beta ; alpha_prime ; beta_prime ; alpha_double_prime ; beta_double_prime];
    x_detrended_mod(i) = inv(I + lmbd^2*(D2'*D2) + gm^2*(G'*G))*(x(i)' + gm^2*G'*gamma);
end
% w = LPFilter(w, 80/fs);

tt = (0:length(x)-1)/fs;

% % figure
% % hold on
% % plot(t, x);
% % plot(t(II), x(II), 'r.');
% % grid
%
% % figure
% % hold on
% % plot(z);
% % plot(II, z(II), 'ro');
% % grid
%
figure
hold on
plot(tt, s, 'k', 'linewidth', 2);
plot(tt, x);
plot(tt, y, 'r');
plot(tt, z, 'g');
plot(tt, w, 'm');
% plot(tt, w, 'm.');
grid

figure
hold on
plot(x - y, 'r');
plot(x - z, 'g');
plot(x - w, 'm');
grid

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2))
snr1 = 10*log10(sum(s.^2)/sum((s - y).^2))
snr2 = 10*log10(sum(s.^2)/sum((s - z).^2))
snr3 = 10*log10(sum(s.^2)/sum((s - w).^2))
snr_detrended = 10*log10(sum(s.^2)/sum((s - x_detrended).^2))
snr_detrended_mod = 10*log10(sum(s.^2)/sum((s - x_detrended_mod).^2))

peaks_s = PeakDetection(s, f0/fs);
I_s = find(peaks_s);
peaks_x = PeakDetection(x, f0/fs);
I_x = find(peaks_x);
peaks_w = PeakDetection(w, f0/fs);
I_w = find(peaks_w);

figure
subplot(211);
plot(tt, s, 'k');
hold on
plot(tt(I_x), s(I_x), 'bo');
plot(tt(I_w), s(I_w), 'ro');
grid

subplot(212);
HR_s = 60*fs./diff(I_s);
HR_x = 60*fs./diff(I_x);
HR_w = 60*fs./diff(I_w);
plot(tt(I_s(2:end)), HR_s, 'ko');
hold on
plot(tt(I_x(2:end)), HR_x, 'b');
plot(tt(I_w(2:end)), HR_w, 'r');
grid

figure
hold on;
plot(x);
plot(x_detrended, 'r', 'linewidth', 2);
plot(x_detrended_mod, 'g', 'linewidth', 2);
grid