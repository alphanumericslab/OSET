function [icasig, W]=milca(x,P,kneig,algo,harm,finetune)

% Mutual Information based Least Dependent Component Analysis
%References:
% H. Stogbauer, A. Kraskov, S. A. Astakhov, P.Grassberger, 
% Phys. Rev. E 70 (6)  066123, 2004
% A. Kraskov, H. Stogbauer, P. Grassberger,  
% Phys. Rev. E 69 (6) 066138, 2004 

% Output: most independent components under linear transformation and the transformation matrix
% x....input data MxN   M...channelnummer  N...sampling points  m<<n
% P - the number of components to be retained in PCA (estimated number of pure sources)
% kneig... k nearest neighbor for MI algorithm
% algo ... version of the MI estimator: 1=cubic  2=rectangular
% harm ... fit the MI vs. angle curve with # harmonics
% finetune ... Number 2^# of angles between 0 and pi/2 (default is 7)

[M,N] = size(x);

%default-values
if ~exist('P'), P=M; end
if ~exist('kneig'), kneig=12; end
if ~exist('algo'), algo=2; end
if ~exist('harm'), harm=1; end
if ~exist('finetune')
    Nb=128;
else
    Nb=2^finetune;
end

% PCA - prewhitening + dimension reduction

covX = cov(x'); 
[E_,D_] = eig(covX);
Whitening_mat = inv(sqrtm(D_((M-P+1):M,(M-P+1):M)))*E_(:,(M-P+1):M)'; 
x=Whitening_mat*x;

% save data for external Programm

zwsp=x';
save zwspmilca.txt zwsp -ASCII


% execute C Programm
[a unout]=system(['milca zwspmilca.txt ',num2str(P),' ',num2str(N),' ',num2str(kneig),' ',num2str(algo-1),' ',num2str(Nb),' ',num2str(harm)]);
Rall=str2num(unout);

%output
icasig=Rall*x;
W=Rall*Whitening_mat;


