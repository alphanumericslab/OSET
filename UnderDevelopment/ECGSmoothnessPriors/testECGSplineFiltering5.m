clear
close all
clc

load ECG2
s = data(4, 1:end);

smoothingwlen1 = 3;
smoothingwlen2 = 3;
smoothingwlen3 = 3;

wlen1 = 15e-3;
wlen2 = 15e-3;
wlen3 = 15e-3;
p1 = 4; % spline order
p2 = 4; % spline order
p3 = 4; % spline order
snr = 25; % dB

s = s - LPFilter(s, 0.1/fs);

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
    t0 = [1 ; zeros(p1, 1)];
    tn = zeros(p1+1, 1);
    t = zeros(p1+1, length(i));
    for k = 1 : p1+1,
        tn(k) = ti(end)^(k-1);
        t(k, :) = ti.^(k-1);
    end
    
    Phi = [t0 tn];
    T = t*t';
    disp(['size(T) = ', num2str(length(T)), ', rank(T) = ', num2str(rank(T))]);
    %     S = Phi'*inv(T)*Phi;
    S = Phi'/T*Phi;
    u = t*x(i)';
    
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
for n = 1 : N - 1,
    i = II(n) : II(n+1); % I used this line to preserve the continuity of the curves
    ti = (0:length(i)-1)/fs;
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
    %     S = Phi'*inv(T)*Phi;
    S = Phi'/T*Phi;
    u = t*x(i)'; % the second stage is still applied to the original data; but we get the boundary conditions from the first stage
    
    %     gamma = [x(i(1)) x(i(end))]';
    % this is an intutive tweak, to force the curve pass an average point
    % rather than the exact signal points at the curve ends:
    alpha = mean(y(max(i(1)-smoothingwlen2,1) : min(i(1)+smoothingwlen2, length(y))));
    beta = mean(y(max(i(end)-smoothingwlen2,1) : min(i(end)+smoothingwlen2, length(y))));
    alpha_prime = mean(fs*diff(y(max(i(1)-smoothingwlen2,1) : min(i(1)+smoothingwlen2, length(y)))));
    beta_prime = mean(fs*diff(y(max(i(end)-smoothingwlen2,1) : min(i(end)+smoothingwlen2, length(y)))));
    gamma = [alpha ; beta ; alpha_prime ; beta_prime];
    
    %     lambda = inv(S)*(Phi'*inv(T)*u - gamma);
    %     aa = inv(T)*(u - Phi*lambda);
    lambda = S\(Phi'/T*u - gamma);
    aa = T\(u - Phi*lambda);
    z(i) = aa'*t;
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
N = length(II);
for n = 1 : N - 1,
    i = II(n) : II(n+1); % I used this line to preserve the continuity of the curves
    ti = (0:length(i)-1)/fs;
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
    S = Phi'/T*Phi;
    u = t*x(i)'; % the third stage is still applied to the original data; but we get the boundary conditions from the first stage
    
    %     gamma = [x(i(1)) x(i(end))]';
    % this is an intutive tweak, to force the curve pass an average point
    % rather than the exact signal points at the curve ends:
    alpha = mean(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))));
    beta = mean(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))));
    alpha_prime = mean(fs*diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z)))));
    beta_prime = mean(fs*diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z)))));
    alpha_double_prime = mean(fs^2*diff(z(max(i(1)-smoothingwlen3,1) : min(i(1)+smoothingwlen3, length(z))), 2));
    beta_double_prime = mean(fs^2*diff(z(max(i(end)-smoothingwlen3,1) : min(i(end)+smoothingwlen3, length(z))), 2));
    gamma = [alpha ; beta ; alpha_prime ; beta_prime ; alpha_double_prime ; beta_double_prime];
    
    %     lambda = inv(S)*(Phi'*inv(T)*u - gamma);
    %     aa = inv(T)*(u - Phi*lambda);
    lambda = S\(Phi'/T*u - gamma);
    aa = T\(u - Phi*lambda);
    w(i) = aa'*t;
end
% w = LPFilter(w, 80/fs);

tt = (0:length(x)-1)/fs;

% figure
% hold on
% plot(t, x);
% plot(t(II), x(II), 'r.');
% grid

% figure
% hold on
% plot(z);
% plot(II, z(II), 'ro');
% grid

figure
hold on
plot(tt, x);
plot(tt, y, 'r');
plot(tt, z, 'k');
plot(tt, w, 'm');
% plot(tt, w, 'm.');
grid

figure
hold on
plot(x - y, 'r');
plot(x - z, 'k');
plot(x - w, 'm');
grid

snr0 = 10*log10(sum(s.^2)/sum((s - x).^2))
snr1 = 10*log10(sum(s.^2)/sum((s - y).^2))
snr2 = 10*log10(sum(s.^2)/sum((s - z).^2))
snr3 = 10*log10(sum(s.^2)/sum((s - w).^2))

