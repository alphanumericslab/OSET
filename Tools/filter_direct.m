function [y yf] = filter_direct(b, a, x, yi)

N = length(x);
y = zeros(1, N);
bb = b(end : -1 : 1);
aa = a(end : -1 : 1);
for i = 1 : N,
    s = 0;
    for m = 1 : length(b)
        s = s + b(m)*x(i - m);
    end
    for n = 1 : length(aa)
        ind = 
        s = s - aa(m)*yy;
    end
    
end

forward = [l ; zeros(length(indexes),1)];
for m = 1 : length(indexes),
    for p = 2 : length(p_causal),
        forward(m + length(l)) = forward(m + length(l)) - p_causal(p)*forward(m + length(l) - (p-1));
    end
    forward(m + length(l)) = (forward(m + length(l)) + sum(p_causal)*xx(m)) / p_causal(1);
end
forward = forward(length(l)+1 : end);
