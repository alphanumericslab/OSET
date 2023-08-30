function [icasig, W]=milcadelay(x,kneig,algo,emb_dim,emb_tau,harm,finetune,PCAflag)

% Mutual information based least dependent component analysis including
% additional the time structure of the input
% Output: most independent components under linear transformation and the transformation matrix
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm
% algo ... 1=cubic  2=rectangular
% emb_dim... embxding dimension (default is 1, no embedding)
% emb_tau... time-delay, only relevant when emb_dim>1 (default is 1)
% harm ... fit the MI vs. angle curve with # harmonics
% finetune ... Number 2^# of angles between 0 and pi/2 (default is 7)
% PCAflag ... 0=without whitening  else whitening  (default with whitening)

%default-values
if ~exist('kneig'), kneig=18; end
if ~exist('algo'), algo=2; end
if ~exist('emb_dim'), emb_dim=2; end
if ~exist('emb_tau'), emb_tau=1; end
if ~exist('harm'), harm=1; end
if ~exist('finetune')
    Nb=128;
else
    Nb=2^finetune;
end
if ~exist('PCAflag'), PCAflag=1; end

[Nd,N]=size(x);
if Nd>N
    x=x';
    [Nd,N]=size(x);
end

%removing mean and normalization
x=x-mean(x')'*ones(1,N); 
x=x./(std(x')'*ones(1,N));

%whitening
if PCAflag
    [u,s,v]=svd(cov(x'));
    V=v*s^(-.5)*u';
    x=V*x;
end

% Augmentation see
% Erik Learnx-Miller and John Fisher, (2003) "ICA Using Spacings Estimates of Entropy." Journal of Machine Learning Research
% x=randn(Nd,N*30)*0.175+repmat(x,[1,30]);

% save data for external Programm
zwsp=x';
save zwspmilcadelay.txt zwsp -ASCII

% execute C Programm
[a unout]=unix(['milcadelay zwspmilcadelay.txt ',num2str(Nd),' ',num2str(N),' ',num2str(kneig),' ',num2str(emb_tau),' ',num2str(emb_dim),' ',num2str(algo-1),' ',num2str(Nb),' ',num2str(harm)]);
Rall=str2num(unout);

%output
icasig=Rall*x;
if PCAflag
    W=Rall*V;
else
    W=Rall;
end



