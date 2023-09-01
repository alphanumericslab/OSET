function MIClustering(x,kneig,algo)

% hierarchical clustering using the grouping property of mutual information
% output: dendrogram 
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm
% algo ... 1=cubic  2=rectangular

%default-values
if ~exist('kneig'), kneig=6; end
if ~exist('algo'), algo=2; end

% check input data if format is correct
[Nd,N]=size(x);
if Nd>N
    x=x';
    [Nd,N]=size(x);
end

zwsp=x';
save zwspMIC.txt zwsp -ASCII

% execute C Programm
[a b]=unix(['MIClustering zwspMIC.txt ',num2str(Nd),' ',num2str(N),' ',num2str(kneig),' ',num2str(algo)]);
clutr=str2num(b);

% Plot Output
figure;
dendrogram(clutr(:,[1,2,5]));
ylabel('mutual information');
title('Hierarchical clustering using the grouping property of mutual information');


