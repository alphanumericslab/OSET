function [x_filtered1, x_filtered2] = ECGSmoothingPriors(x, DiffOrder, wlen, lambda, guardlen_l, guardlen_u, adapt, withwindow)

x = x(:);
knots1 = 1: round(wlen) : length(x);
% add the first and last samples as knots to avoid end-point discontinuities
if(knots1(1) > 1)
    knots1 = [1 knots1];
end
if(knots1(end) < length(x))
    knots1 = [knots1 length(x)];
end

% first round
x_filtered1 = zeros(size(x));
for n = 1 : length(knots1) - 1
    i = knots1(n) : knots1(n+1);
    ll = length(i);
    I = speye(ll);
    Dd = spdiags(ones(ll-2, 1)*diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, ll-DiffOrder, ll);
    %     R = qr(I + lambda^2*(Dd'*Dd));
    %     x_filtered1(i) = R\(R'\x(i));
    if(withwindow)
        W = diag(hamming(length(I)));
    else
        W = I;
    end
    x_filtered1(i) = (W + lambda^2*(Dd'*Dd))\(W*x(i));
end

% dd = full(Dd);
% dd2 = dd'*dd;
% dd2 = inv(eye(length(dd2)) + lambda^2*dd2);
% dd2(1:10,1:10);

% second round
knots2 = round(wlen/2): round(wlen) : length(x);
% II = II - round(wlen*fs/2);II(II < 1) = [];
if(knots2(1) > 1)
    knots2 = [1 knots2];
end
if(knots2(end) < length(x))
    knots2 = [knots2 length(x)];
end

x_filtered2 = zeros(size(x));
for n = 1 : length(knots1) - 1
    i = knots2(n) : knots2(n+1); % I used this line to preserve the continuity of the curves
    
    ll = length(i);
    I = speye(ll);
    Dd = spdiags(ones(ll -2 + guardlen_l + guardlen_u, 1)*diff([zeros(1, DiffOrder) 1 zeros(1, DiffOrder)], DiffOrder), 0:DiffOrder, ll + guardlen_l + guardlen_u - DiffOrder, ll + guardlen_l + guardlen_u);
    D_l = Dd(:, 1:guardlen_l);
    D = Dd(:, guardlen_l+1:end-guardlen_u);
    D_u = Dd(:, end-guardlen_u+1:end);
    indxl = i(1)-guardlen_l : i(1)-1;
    indxl(indxl < 1) = 1;
    
    indxu = i(end) + 1 : i(end) + guardlen_u;
    indxu(indxu > length(x_filtered1)) = length(x_filtered1);
    
    l = x_filtered1(indxl);
    u = x_filtered1(indxu);
    if(adapt)
        lambda_adaptive = (lambda/(sqrt(length(x)/length(x(i)))*norm(x(i))/norm(x))); % adaptive version
    else
        lambda_adaptive = lambda; % fixed version
    end
    
    if(withwindow)
        W = diag(hamming(length(I)));
    else
        W = I;
    end
    %     R = qr(I + lambda_adaptive^2*(D'*D));
    %     x_filtered2(i) = R\(R'\(x(i) - lambda_adaptive^2*(D'*D_l(:, 1:length(l)))*l - lambda_adaptive^2*(D'*D_u(:, 1:length(u)))*u));
    x_filtered2(i) = (W + lambda_adaptive^2*(D'*D))\(W*x(i) - lambda_adaptive^2*(D'*D_l(:, 1:length(l)))*l - lambda_adaptive^2*(D'*D_u(:, 1:length(u)))*u);
end

x_filtered1 = x_filtered1(:)';
x_filtered2 = x_filtered2(:)';
