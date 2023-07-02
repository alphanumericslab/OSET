function y = RemoveEEGSeizureNoise(x, fsx, probth, wlen, keep)

% NSCA
lx = round(wlen*fsx);
Ex = sqrt(filter(ones(1, lx), lx, x.^2, [], 2));
eex = mean(Ex, 1);
thx = quantile(eex, probth);
Ix = find(eex >= thx);
Jx = 1:length(Ex); %find(eex < thx);%
[s , ~, A] = NSCA(x, Ix, Jx);
y = A(:, end-keep+1:end)*s(end-keep+1:end, :);
