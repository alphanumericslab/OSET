function out=amari(C,A)

[b,a]=size(C);

dummy=pinv(A)*C;
dummy=sum(ntu(abs(dummy)))-1;

dummy2=pinv(C)*A;
dummy2=sum(ntu(abs(dummy2)))-1;

out=(sum(dummy)+sum(dummy2))/(2*a*(a-1));


function CN=ntu(C)
[m n]=size(C);
for t=1:n,
 CN(:,t)=C(:,t)./max(abs(C(:,t)));
end
