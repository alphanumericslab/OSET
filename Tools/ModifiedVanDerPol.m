function xout = ModifiedVanDerPol(x)
% An implementation of the modified van Der Pol equation

w0 = 1;
alpha = .75;
a = .1;
b = .5;
c = .9;

% % % xout = [0 1 ; -w0^2 -2*alpha]*x(:);
% % % xout = [0 1 ; -w0^2 -2*(x(1)^2-1)*(x(1)^2-3)*(x(1)^2-5)]*x(:);
% % % xout = [0 1 ; -w0^2 -2*alpha*(x(1)^2-1)]*x(:);
% % % xout = [0 1 ; -w0^2 2*alpha*(abs(x(1))-1)*(abs(x(1))-1.5)*(abs(x(1))-2)]*x(:);
% % % xout = [0 1 ; -w0^2 -(x(1)^2 - 1)*(x(1)^2 - 1.5)*(x(1)^2 - 2)]*x(:);
% xout = [0 1 ; -w0^2 -alpha*(x(1)^2 - a^2)*(x(1)^2 - b^2)*(x(1)^2 - c^2)]*x(:);
xout = [0 1 ; -w0^2 -alpha*(abs(x(1)) - a)*(abs(x(1)) - b)*(abs(x(1)) - c)]*x(:);

% % % zeta = .9;
% % % w0 = (x(1)^2 - 1)*(x(1)^2 - 1.5)*(x(1)^2 - 2);
% % % xout = [0 1 ;  -w0^2 -2*w0*zeta]*x(:);
