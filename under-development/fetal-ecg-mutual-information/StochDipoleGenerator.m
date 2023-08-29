function [DIP teta]= StochDipoleGenerator(N,fs,f,alphai,bi,tetai,fSD,alphaiSD,biSD,tetaiSD,teta0);

w = 2*pi*f;
dt = 1/fs;

X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);
teta = zeros(1,N);
teta(1,1:end) = [(teta0):w*dt:((w*dt)*(N-1) + teta0)] + 2*pi*fSD*dt*randn(1,N);

for j = 1:length(alphai.x),
    dtetaix = rem(teta - tetai.x(j) + pi + randn(1,N)*tetaiSD.x(j),2*pi)-pi;
    X = X + (alphai.x(j)+randn(1,N)*alphaiSD.x(j)) .* exp(-dtetaix .^2 ./ (2*(bi.x(j)+randn(1,N)*biSD.x(j)) .^ 2));   % x state variable
end

for j = 1:length(alphai.y),
    dtetaiy = rem(teta - tetai.y(j) + pi + randn(1,N)*tetaiSD.y(j),2*pi)-pi;
    Y = Y + (alphai.y(j)+randn(1,N)*alphaiSD.y(j)) .* exp(-dtetaiy .^2 ./ (2*(bi.y(j)+randn(1,N)*biSD.y(j)) .^ 2));   % y state variable
end

for j = 1:length(alphai.z),
    dtetaiz = rem(teta - tetai.z(j) + pi + randn(1,N)*tetaiSD.z(j),2*pi)-pi;
    Z = Z + (alphai.z(j)+randn(1,N)*alphaiSD.z(j)) .* exp(-dtetaiz .^2 ./ (2*(bi.z(j)+randn(1,N)*biSD.z(j)) .^ 2));   % z state variable
end

DIP.x = X;
DIP.y = Y;
DIP.z = Z;
teta = rem(teta+ pi,2*pi)-pi;
