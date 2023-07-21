%//////////////////////////////////////////////////////////////////////////
function C = Cost_gamma(x, b, s, sigma2, gamma)
% Cost_gamma has been deprecated and integrated with OptimalSmoothnessFactor.
% n = min([length(x), length(b), length(s)]);
% x = x(1:n);
% b = b(1:n);
% s = s(1:n);
% C = sum(((s.^2.*x + s.*b)./(gamma + s.^2)).^2) - sigma2;
% end
error('Cost_gamma has been deprecated and integrated into OptimalSmoothnessFactor.');
