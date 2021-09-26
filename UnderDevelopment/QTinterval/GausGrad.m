function gr=GausGrad(t,p)
gr=zeros(length(p),length(t)); p=p(:);
for i=1:3:length(p)
    gr(i,:)=GausVal(t,[1 p(i+1) p(i+2)]);
    gr(i+1,:)=(((t-p(i+2)).^2)/(p(i+1)^3)).*GausVal(t,[p(i) p(i+1) p(i+2)]); 
    gr(i+2,:)=((t-p(i+2))/(p(i+1)^2)).*GausVal(t,[p(i) p(i+1) p(i+2)]);
end
