function [y,W,A] = SCAmodified(x,f1,f2) % spectral component analysis
%
% f1,f2: Normalized frequencies

xx = BPFilter(x, f1, f2);
% xx = BPFilter(xx, f1, f2);

% h = fdesign.bandpass('N,F3dB1,F3dB2,Ast1,Ap,Ast2', 120, f1, f2, 50, 0.5, 50);
% Hd = design(h);% fvtool(Hd);
% hdf1 = convert(Hd,'df1');
% xx = filter(hdf1,x,2);

B = xx*xx';
C = x*x';

C = (C+C')/2;
B = (B+B')/2;

[V,D] = eig(B,C,'chol');

d = diag(D);
[~,II] = sort(d);
II = II(end:-1:1);

W = V(:,II)';
A = pinv(W);

y = real(W*x);

