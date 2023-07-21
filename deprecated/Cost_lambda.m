%//////////////////////////////////////////////////////////////////////////
function C = Cost_lambda(x, b, s, epsilon2, lambda)
% Cost_lambda has been deprecated and integrated into OptimalSmoothnessFactor.
% function C = Cost_lambda(x, b, s, epsilon2, lambda)
% n = min([length(x), length(b), length(s)]);
% x = x(1:n);
% b = b(1:n);
% s = s(1:n);
% C = sum(((s.*x + b)./(1 + lambda.*s.^2)).^2) - epsilon2;
% end
error('Cost_lambda has been deprecated and integrated with OptimalSmoothnessFactor.');

