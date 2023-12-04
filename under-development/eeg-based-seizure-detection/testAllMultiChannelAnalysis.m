close all;
clear;
clc;

% parameters
subject1 = 'Dog_4';
subject2 = 'Dog_4';
tri = 28;
number1 = tri;
number2 = tri;
path1 = ['E:\Sameni\Projects\Seizure\' subject1 '\'];
path2 = ['E:\Sameni\Projects\Seizure\' subject2 '\'];
mode1 = 'interictal';
mode2 = 'preictal';
SpectralWinLen = 5.0; % in seconds
SpectralOverlapPercentage = 90; % [0 99]
ffs = 120;
f0 = 60; % powerline frequency Hz
Q = 60; % notch filter Q-factor
nfft = 256;
wlen = 2.0; % energy calculation window length in seconds
% wlen2 = 3.0; % second energy calculation window length in seconds
probthx = 0.8; % quantile calculation threshold
probthy = 0.8; % quantile calculation threshold
qtl = 0.75; % the samples above this quantile are excluded from NSCA calculations
wdeclevels = 10;
% keep = 3;
sourcestodenoise = 8;
drop = 500; % the number of samples to exclude from the head and tail of the data during processing
medianwin = 0.025; % median calculation window length in seconds
median_sigma_factor = 400.0;% used in the median filter threshold
wlow = 1; % used in WNSCA
whigh = 2.5; % used in WNSCA
sourcenum = 12;
epochlen = 60.0; % epoch length in seconds
epochoverlap = 55.0; % overlap between epochs in seconds

% load signals
[x0 fsx0] = LoadSeizureEEGFull(path1, subject1, mode1, number1);
[y0 fsy0] = LoadSeizureEEGFull(path2, subject2, mode2, number2);

% preprocessing: resample the data to 120Hz with a polyphase bandpass filter passing [1Hz-55Hz]
load BandpassFilter1Hzto55HzAt400Hz h % load the bandpass filter used for resampling the input data
if(fsx0 < 500)
    x = resample(x0', ffs, round(fsx0), h*round(fsx0)/ffs)';
else
    midfs = 400;
    % resample to intermediate frequency
    x_ = resample(x0', midfs, round(fsx0))';
    
    % IIR notch
    Wo = f0/(midfs/2);  BW = Wo/Q;
    [b,a] = iirnotch(Wo, BW);
    x_ = filter(b, a, x_, [], 2);
    
    % resample using bandpass filter
    x = resample(x_', ffs, midfs, h*midfs/ffs)';
    
    clear x_ midfs b a Wo BW
end
fsx = ffs;
if(fsy0 < 500)
    y = resample(y0', ffs, round(fsy0), h*round(fsy0)/ffs)';
else
    midfs = 400;
    % resample to intermediate frequency
    y_ = resample(y0', midfs, round(fsy0))';
    
    % IIR notch
    Wo = f0/(midfs/2);  BW = Wo/Q;
    [b,a] = iirnotch(Wo, BW);
    y_ = filter(b, a, y_, [], 2);
    
    % resample using bandpass filter
    y = resample(y_', ffs, midfs, h*midfs/ffs)';
    
    clear y_ midfs b a Wo BW
end
fsy = ffs;

% % just for test
% x = LPFilter(x, 10.0/fsx);
% y = LPFilter(y, 10.0/fsx);

% drop transient samples
x = x(:,drop:end-drop);
y = y(:,drop:end-drop);

% normalize
x = (x - mean(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
y = (y - mean(y, 2)*ones(1, size(y,2)))./(std(y, [], 2)*ones(1, size(y,2)));
% x = (x - median(x, 2)*ones(1, size(x,2)))./(std(x, [], 2)*ones(1, size(x,2)));
% y = (y - median(y, 2)*ones(1, size(y,2)))./(std(y, [], 2)*ones(1, size(y,2)));

bins = 256;
Hx = zeros(size(x,1), bins);
Hy = zeros(size(y,1), bins);
Hvalsx = zeros(size(x,1), bins);
Hvalsy = zeros(size(y,1), bins);
for i = 1:size(x, 1)
    [Hx(i, :) Hvalsx(i, :)] = hist(x(i, :), bins);
end
for i = 1:size(y, 1)
    [Hy(i, :) Hvalsy(i, :)] = hist(y(i, :), bins);
end

cdfx = cumsum(Hx, 2); cdfx = cdfx./(cdfx(:,end)*ones(1, size(cdfx,2)));
cdfy = cumsum(Hy, 2); cdfy = cdfy./(cdfy(:,end)*ones(1, size(cdfy,2)));

xx = zeros(size(x));
for i = 1:size(x, 1)
    for j = 1:size(x, 2)
        I = find(Hvalsx(i, :) >= x(i, j), 1, 'first');
        if(isempty(I))
            I = size(cdfx,2);
        end
        xx(i, j) = cdfx(i, I) - 0.5;
    end
end
yy = zeros(size(y));
for i = 1:size(y, 1)
    for j = 1:size(y, 2)
        I = find(Hvalsy(i, :) >= y(i, j), 1, 'first');
        if(isempty(I))
            I = size(cdfy,2);
        end
        yy(i, j) = cdfy(i, I) - 0.5;
    end
end

xx = x;
yy = y;

% PCA
Cx = cov(xx');
[Vx Dx] = eig(Cx);
dx = diag(Dx);
[~,Ix] = sort(dx, 1, 'descend');
wxPCA = Vx(:, Ix(1:sourcenum))';
% qx = wxPCA*xx;
qx = wxPCA*x;
dx = dx(Ix)/dx(Ix(1));

Cy = cov(yy');
[Vy Dy] = eig(Cy);
dy = diag(Dy);
[~,Iy] = sort(dy, 1, 'descend');
wyPCA = Vy(:, Iy(1:sourcenum))';
% qy = wyPCA*yy;
qy = wyPCA*y;
dy = dy(Iy)/dy(Iy(1));

% JADE
wxJADE = jadeR(xx, sourcenum);
% sx = wxJADE*xx;
sx = wxJADE*x;
wyJADE = jadeR(yy, sourcenum);
% sy = wyJADE*yy;
sy = wyJADE*y;

% SCA
[~, wxSCA] = SCA2(xx, 5.0/fsx, 8.0/fsx, 5);
[~, wySCA] = SCA2(yy, 5.0/fsy, 8.0/fsy, 5);
px = wxSCA*x;
py = wySCA*y;

% NSCA
lx = round(wlen*fsx);
ly = round(wlen*fsy);
Ex = sqrt(filter(ones(1, lx), lx, xx.^2, [], 2));
Ey = sqrt(filter(ones(1, ly), ly, yy.^2, [], 2));
eex = mean(Ex, 1);
eey = mean(Ey, 1);
thx = quantile(eex, probthx);
thy = quantile(eey, probthy);
Ix = find(eex >= thx);
Jx = 1:length(Ex);%find(eex < thx);%
Iy = find(eey >= thy);
Jy = 1:length(Ey);%find(eey < thy);%
[~ , wxNSCA] = NSCA(xx, Ix, Jx);
[~ , wyNSCA] = NSCA(yy, Iy, Jy);
fx = wxNSCA*x;
fy = wyNSCA*y;

% WNSCA
% Ix = 1:length(xx);
% Jx = 1:length(xx);
wxI = linspace(wlow, whigh, length(Ix));
% wxI = logspace(0, 1, length(Ix));
wxJ = linspace(1, 1, length(Jx));
[~ , wxWNSCA2] = WNSCA2(xx, Ix, wxI, Jx);
% Iy = 1:length(yy);
% Jy = 1:length(yy);
wyI = linspace(wlow, whigh, length(Iy));
% wyI = logspace(0, 1, length(Iy));
wyJ = linspace(1, 1, length(Jy));
[~ , wyWNSCA2] = WNSCA2(yy, Iy, wyI, Jy);
gx = wxWNSCA2*x;
gy = wyWNSCA2*y;

% PlotECG(x, 8, 'c', fsx);
% PlotECG(xx, 8, 'k', fsx);
PlotECG(qx(1:sourcenum, :), 4, 'b', fsx, 'PCA');
PlotECG(sx(1:sourcenum, :), 4, 'b', fsx, 'JADE');
PlotECG(px(1:sourcenum, :), 4, 'b', fsx, 'SCA');
PlotECG(fx(1:sourcenum, :), 4, 'b', fsx, 'NSCA');
PlotECG(gx(1:sourcenum, :), 4, 'b', fsx, 'WNSCA');

% PlotECG(y, 8, 'g', fsy);
% PlotECG(yy, 8, 'm', fsy);
PlotECG(qy(1:sourcenum, :), 4, 'r', fsy, 'PCA');
PlotECG(sy(1:sourcenum, :), 4, 'r', fsy, 'JADE');
PlotECG(py(1:sourcenum, :), 4, 'r', fsy, 'SCA');
PlotECG(fy(1:sourcenum, :), 4, 'r', fsy, 'NSCA');
PlotECG(gy(1:sourcenum, :), 4, 'r', fsy, 'WNSCA');

