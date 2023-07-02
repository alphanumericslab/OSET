function [y, B, conv] = EASI(x, nsou, lambda, nlintype)
% A simple demo program for th EASI algorithm:
% an adaptive blind source separation algorithm
% Adopted from:
% Ref :
%   @inproceedings{PFSEusipco,
%   title = "Adaptive source separation without prewhitening",
%   author = "Beate Laheld and Jean-Fran\c{c}ois Cardoso",
%   booktitle = "{Proc. EUSIPCO'94}",
%   pages = "183--186",
%   address = {Edinburgh},
%   month = sep,
%   year = "1994" }
% NOTE: This code has been adapted from a demo file obtained from:
% http://bsp.teithe.gr/members/downloads/easidemo/easi_demo.m

ncap = size(x, 1);
T = size(x, 2);
idsou	= eye(nsou);

B = randn(nsou,ncap);
y = zeros(nsou, T);
conv = zeros(1, T);
for t = 1:T
    y(:, t)	= B * x(:, t);
    y2	= y(:, t) * y(:, t)';
    switch nlintype
        case 1, g	= diag(y2) .* y(:, t);
        case 2, g	=   y(:, t) .* diag(0.1*y2).*diag(y2);
        case 3, g	=   y(:, t) .* sqrt(diag(y2));
        case 4, g	=   y(:, t) .* log(diag(y2));
        case 5, g	= - y(:, t) .* (diag(y2)<0.9);
        case 6, g	=   y(:, t) ./ log(diag(y2));
        case 7, g	= - y(:, t) ./ sqrt(diag(y2));
        case 8, g	= - y(:, t) ./ diag(y2);
    end
    gy	= g * y(:, t)';
    G	= (y2 - idsou)/(1+lambda*trace(y2)) + (gy - gy')/(1+lambda*abs(g'*y(:, t)));
    B	= B - lambda*G*B;
    conv(t) = norm(G);
end