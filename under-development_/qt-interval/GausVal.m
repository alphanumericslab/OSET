function g=GausVal(t,p)
t=t(:);
if isvector(p)
    p=reshape(p,3,[]);
end
g=sum(p(1,:).*exp(-((repmat(t,1,size(p,2))-p(3,:)).^2)./(2*p(2,:).^2)),2);
end
