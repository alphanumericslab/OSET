function y = BPFilterComplex(x,w0,N)
% A bandpass CIC filter
% 
% Reza Sameni
% December 2008

L1 = size(x,1);
L2 = size(x,2);

z = x .* (ones(L1,1)*exp(-1j*2*pi*w0*(0:L2-1)));
u = filter(ones(1,N), N, z')';
y = u .* (ones(L1,1)*exp(1j*2*pi*w0*(0:L2-1)));