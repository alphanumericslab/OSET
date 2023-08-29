function R = RotateNDim(N,axes,tetas)

if (size(axes,1)~=2 || size(axes,2)~=length(tetas))
    error('Dimension mismatch');
end
if N*(N-1)/2<length(tetas)
    error('Too many rotation angles');
end

L = length(tetas);
R = eye(N);
for i = 1:L
    R(axes(1,i),axes(1,i)) = cos(tetas(i));
    R(axes(2,i),axes(2,i)) = cos(tetas(i));
    R(axes(1,i),axes(2,i)) = sin(tetas(i));
    R(axes(2,i),axes(1,i)) = -sin(tetas(i));
end