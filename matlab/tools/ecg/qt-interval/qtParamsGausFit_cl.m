function [GaussParams, varargout]=qtParamsGausFit_cl(ecg, fs, varargin)
% 'corrected levels' version.
% estimate the Q-T interval and some other Q-T parameters by fitting 
% Gaussian function models.
% 
% Syntax:
% GaussParams=qtParamsGausFit(ecg, fs)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi, lb, ub)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi, lb, ub, beta)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi, lb, ub, beta, options)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi, lb, ub, beta, PrMu, PrCov, varNoise)
% GaussParams=qtParamsGausFit(ecg, fs, rpeaks, p0, soi, lb, ub, beta, PrMu, PrCov, varNoise, options)
% [GaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit(ecg, fs, ...)
% 
% 
% 
% Inputs:
% ecg: matrix of multichannel ECG data. Each column corresponds to each channel.
% fs: sampling frequency.
% rpeaks (or hr): if vector, it is index vector of R peak positions with the
%   same lenth of ecg columns. rpeaks has 1 or true values at R peaks and 0
%   or flase values elsewhere. 
%   if scalar, it is the average heart rate in Hz, used in R peak detection.
% p0: a structure in which p0.q and p0.t are matrices of parameters' initial 
%   values. Their each column corresponds to each channels, and is like
%   [a1 b1 c1 a2 b2 c2 ...]' where ai, bi and ci are respectively the i-th 
%   Gaussians' amplitudes, widths (sec), and centers (sec). R peak is the reference time
%   for c.
% soi: segments-of-interest start and end time (sec), in a structure format in 
%   which soi.q = [q_start  q_end] and soi.t = [t_start  t_end]. R peak is the 
%   reference time (=0) in each beat. rpeaks must be provided in inputs when you
%   use soi. Also, soi.q and soi.t can be a Nx2 matrix where N is the
%   number of beats (length of rpeaks), in which each row contains the soi
%   of each beat.
% lb: lower bound for Gaussians parameters with the same format as p0.
% ub: upper bound for Gaussians parameters with the same format as p0.
% beta: (scalar) constant coefficient in (wave width)=2*beta*(Gaussian width), and 
%   wave offset/onset = Gaussian center) +/- beta*(Gaussian width). In
%   order to define it for Q and T waves, use a structure form like beta.q
%   and beta.t.
% PrMu: mean of prior distribution with the same format as p0.
% PrCov: a structure containing covarinace matrices of prior information. 
%   PrCov.q is covariace matrix of Q wave quassian parameters, and PrCov.t 
%   is for T wave. Each matrix is 3G-by-3G-by-C where G is the 
%   number of Gaussians, and C is the number of channels.
% varNoise: variance of the additive noise
% options: options structure for optimization problem. (Note: It should be in a
%   simple strutre format, NOT the 'optimoptions' format commonly used in matlab
%   optimization toolbox. However, the fields are same as 'optimoptions'.)
% 
% Note: each of rpeaks, p0, soi, lb, ub, beta can be left empty ([]), so the
%   default values will be used
% 
% Output:
% GaussParams: a structure containing q and t fields, for Q and T waves' 
%   estimated parameters of Gaussians, in a 3G-by-B-by-C matrix where G is
%   the number of Gaussians, B is the number of beats, and C in the number 
%   of channels. 
% rPeaks: the R peak positions (sample)
% waveParams: the onset anf offset of the Q and T waves, in waveParams.q
%    and waveParams.t respectively.
% qtInt: Q-T intervals lendth.
% 
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

%% initial values
options=struct; rPeaks=[]; p0.q=[]; p0.t=[]; soi.q=[]; soi.t=[]; 
lb.q=[]; lb.t=[]; ub.q=[]; ub.t=[]; beta.q=[]; beta.t=[]; Bys=false;
hr= 70/60; % Hz, heart rate used in r peak detector
autoInit.q=false; autoInit.t=false;

len = size(ecg,1);

nChannels=size(ecg,2);
%% inputs or default values:
if nargin>=3 && ~isscalar(varargin{1})
    rPeaks=varargin{1};
elseif nargin>=3 && isscalar(varargin{1})
    hr=varargin{1};
end
if length(rPeaks) == len
    rPeaks = find(rPeaks);
end
if isempty(rPeaks)
    [~, rPeaks] = peak_det_amp_threshold(mean(ecg,2),hr/fs, .6, 2, 'MEDIAN');  % R peaks position
end
nBeats=length(rPeaks); 

%%%%%%%%
if nargin>=4
    p0=varargin{2};
end
if isempty(p0)
    p0.q = zeros(3,0);
    p0.t = zeros(3,0);
end
if isempty(p0.q)
    autoInit.q = true;
    p0.q = zeros(3,0);
end
if isempty(p0.t)
    autoInit.t = true;
    p0.t = zeros(3,0);
end
%%%%%%%%
if nargin>=5
    soi=varargin{3};
end
if isempty(soi)
    soi.q = ones(nBeats,1).*[-.05 -.015];
    soi.t = ones(nBeats,1).*[.1 .5];
end
if isempty(soi.q)
    soi.q = ones(nBeats,1).*[-.05 -.015];
end
if isempty(soi.t)
    soi.t = ones(nBeats,1).*[.1 .5];
end

if isvector(soi.q)
    soi.q = ones(nBeats,1).*soi.q(:)';
end
if isvector(soi.t)
    soi.t = ones(nBeats,1).*soi.t(:)';
end


%%%%%%%%
if nargin>=6
    lb=varargin{4};
end
if isempty(lb)
    lb.q=zeros(size(p0.q));
    lb.q(1:3:end,:)=-inf;
    lb.q(2:3:end,:)=.001;
    lb.q(3:3:end,:)=min(soi.q(:,1));
    lb.t=zeros(size(p0.t));
    lb.t(1:3:end,:)=-inf;
    lb.t(2:3:end,:)=.001;
    lb.t(3:3:end,:)=min(soi.t(:,1));
end
if isempty(lb.q)
    lb.q=zeros(size(p0.q));
    lb.q(1:3:end,:)=-inf;
    lb.q(2:3:end,:)=.001;
    lb.q(3:3:end,:)=min(soi.q(:,1));
end
if isempty(lb.t)
    lb.t=zeros(size(p0.t));
    lb.t(1:3:end,:)=-inf;
    lb.t(2:3:end,:)=.001;
    lb.t(3:3:end,:)=min(soi.t(:,1));
end
%%%%%%%
if nargin>=7
    ub=varargin{5};
end
if isempty(ub)
    ub.q=zeros(size(p0.q));
    ub.q(1:3:end,:)=inf;
    ub.q(2:3:end,:)=inf;
    ub.q(3:3:end,:)=max(soi.q(:,2));
    ub.t=zeros(size(p0.t));
    ub.t(1:3:end,:)=inf;
    ub.t(2:3:end,:)=inf;
    ub.t(3:3:end,:)=max(soi.t(:,2));
end
if isempty(ub.q)
    ub.q=zeros(size(p0.q));
    ub.q(1:3:end,:)=inf;
    ub.q(2:3:end,:)=inf;
    ub.q(3:3:end,:)=max(soi.q(:,2));
end
if isempty(ub.t)
    ub.t=zeros(size(p0.t));
    ub.t(1:3:end,:)=inf;
    ub.t(2:3:end,:)=inf;
    ub.t(3:3:end,:)=max(soi.t(:,2));
end
%%%%%%%%%%%
if nargin>=8
    beta=varargin{6};
end
if isempty(beta)
    beta.q=3;
    beta.t=3;
elseif ~isstruct(beta)
    Beta=beta;
    beta=struct;
    beta.q=Beta;
    beta.t=Beta;    
end
if isempty(beta.q)
    beta.q=3;
end
if isempty(beta.t)
    beta.t=3;
end
%%%%%%%
if  nargin>=10 % Bys
    Bys=true;
    PrMu=varargin{7};
    PrCov=varargin{8};
    varNoise=varargin{9};
end
if nargin==9 || nargin==12
    options=varargin{end};
else
    options = struct('SpecifyObjectiveGradient',true);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gaussian fitting
GaussParams.q=nan(size(p0.q,1),nBeats,nChannels);
GaussParams.t=nan(size(p0.t,1),nBeats,nChannels);

for i=1:nBeats
    qOn=floor(soi.q(i,1)*fs);
    qOff=ceil(soi.q(i,2)*fs);
    tOn=floor(soi.t(i,1)*fs);
    tOff=ceil(soi.t(i,2)*fs);
    if rPeaks(i) +  qOn >= 1
        for j=1:nChannels
            %%%% setting initial point 
            if autoInit.q
                [qp, qpi] = min(ecg(limit_indexes(rPeaks(i)+ (qOn:qOff), len) ,j));
                p0.q = [qp; .02;  soi.q(i,1) + qpi/fs];
            end
            %%% dc-level correcion
            isolvl = mean(ecg(limit_indexes(rPeaks(i)+ (qOn - floor(.04*fs):qOn), len),j));
            %%%% gauss fitting
            if Bys
                GaussParams.q(:, i, j)=GausFit((qOn:qOff)/fs, ...
                    ecg(rPeaks(i)+ (qOn:qOff) ,j) - isolvl, p0.q, lb.q, ub.q, PrMu.q(:,j), ...
                    PrCov.q(:,:,j), varNoise ,options);
            else
                GaussParams.q(:, i, j)=GausFit((qOn:qOff)/fs, ...
                    ecg(rPeaks(i)+(qOn:qOff) ,j) - isolvl, p0.q, lb.q, ub.q, options);
            end
        end
    end
    if rPeaks(i) +  tOff <= size(ecg,1)
        for j=1:nChannels
            %%%% setting initial point 
            if autoInit.t
                if abs(max(ecg(rPeaks(i)+ (tOn:tOff) ,j))) > abs(min(ecg(rPeaks(i)+ (tOn:tOff) ,j))) % upward T wave
                    [tp, tpi] = max(ecg(rPeaks(i)+ (tOn:tOff) ,j));
                else % downward T wave
                    [tp, tpi] = min(ecg(rPeaks(i)+ (tOn:tOff) ,j));
                end
                p0.t = [tp; .05; soi.t(i,1) + tpi/fs];
            end
            %%% dc-level correcion
            stlvl = mean(ecg(rPeaks(i)+(floor(.04*fs):floor(.1*fs)),j));
            %%%% gauss fitting
            if Bys
                GaussParams.t(:, i, j)=GausFit((tOn:tOff)/fs, ...
                    ecg(rPeaks(i)+(tOn:tOff),j) - stlvl, p0.t, lb.t, ub.t, PrMu.t(:,j), ...
                    PrCov.t(:,:,j), varNoise ,options);
            else
                GaussParams.t(:, i, j)=GausFit((tOn:tOff)/fs, ...
                    ecg(rPeaks(i)+(tOn:tOff),j) - stlvl, p0.t, lb.t, ub.t, options);
            end
        end
    end
end

varargout{1}=rPeaks;
waveParams.q= [GaussParams.q(3,:,:) - beta.q*GaussParams.q(2,:,:);
     GaussParams.q(3,:,:) + beta.q*GaussParams.q(2,:,:)];
waveParams.t= [GaussParams.t(3,:,:) - beta.t*GaussParams.t(2,:,:);
     GaussParams.t(3,:,:) + beta.t*GaussParams.t(2,:,:)];
varargout{2} = soi;
varargout{3} = waveParams;
t_offset = GaussParams.t(3,:,:) + beta.t*GaussParams.t(2,:,:);
q_onset = GaussParams.q(3,:,:) - beta.q*GaussParams.q(2,:,:);
varargout{4} = t_offset - q_onset;
end

function refined_indexes = limit_indexes(indexes, len)
    first_index = find(indexes > 0 , 1, 'first');
    last_index = find(indexes <= len , 1, 'last');
    refined_indexes = indexes(first_index : last_index);
end
  
