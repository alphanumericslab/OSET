function f = IFMCorrelator(x, f0, BW, tau)
%
% An instantaneous frequency estimator based on lagged correlations
% Reza Sameni
% May 2019

x = x(:)';

T = length(x);
n = (0:T-1);

c = cos(2*pi*f0*n);
s = -sin(2*pi*f0*n);

y_i = x .* c;
y_q = x .* s;

z_i = LPFilter(y_i, BW);
z_q = LPFilter(y_q, BW);

z_i_lagged = [z_i(tau+1:end) zeros(1, tau)];
z_q_lagged = [z_q(tau+1:end) zeros(1, tau)];

A1 = z_i_lagged - z_q;
B1 = z_q_lagged - z_i;
C1 = z_q_lagged + z_q;
D1 = -z_i_lagged + z_i;

A2 = A1.^2;
B2 = B1.^2;
C2 = C1.^2;
D2 = D1.^2;

A3 = LPFilter(A2, BW);
B3 = LPFilter(B2, BW);
C3 = LPFilter(C2, BW);
D3 = LPFilter(D2, BW);

E = B3 - A3;
F = C3 - D3;

f = atan2(E, F)/(2*pi*tau);

