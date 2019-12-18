close all;
clear;
clc;

% parameters
subject1 = 'Dog_3';
subject2 = 'Dog_3';
tri = 9;
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
wlen2 = 3.0; % second energy calculation window length in seconds
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
whigh = 10; % used in WNSCA
sourcenum = 6;
epochlen = 10.0; % epoch length in seconds
epochoverlap = 9.75; % overlap between epochs in seconds

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

% % pre-process
% % x = LPFilter(x, 10.0/fsx); % lowpass filter
% xx = x - LPFilter(x, .25/fsx); % highpass filter
%
% % y = LPFilter(y, 10.0/fsy); % lowpass filter
% yy = y - LPFilter(y, .25/fsy); % highpass filter

% trim the extreme peaks
medx = zeros(size(x));
xx = x;
for i = 1: size(x, 1)%sourcestodenoise,%
    %     xx(i, :) = x(i, :) - wden(x(i, :),  'rigrsure', 's', 'mln', wdeclevels, 'db8');
    %     medx(i, :)  = TrimmedFilter(x(i, :),'trmean', round(medianwin*fsx), ceil(round(medianwin*fsx)/10.0) );
    medx(i, :)  = TrimmedFilter(x(i, :),'median', round(medianwin*fsx) );
    %     ex = abs(x(i, :) - medx(i, :));
    ex = abs(x(i, :));
    medthx = median_sigma_factor*std(ex);
    Ix = (ex > medthx);
    xx(i, Ix) = medx(i, Ix);
end

medy = zeros(size(y));
yy = y;
for i = 1: size(y, 1)%sourcestodenoise,%
    %     yy(i, :) = y(i, :) - wden(y(i, :),  'rigrsure', 's', 'mln', wdeclevels, 'db8');
    %     medy(i, :)  = TrimmedFilter(y(i, :),'trmean', round(medianwin*fsy), ceil(round(medianwin*fsy)/10.0) );
    medy(i, :)  = TrimmedFilter(y(i, :),'median', round(medianwin*fsy) );
    %     ey = abs(y(i, :) - medy(i, :));
    ey = abs(y(i, :));
    medthy = median_sigma_factor*std(ey);
    Iy = (ey > medthy);
    yy(i, Iy) = medy(i, Iy);
end

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
[sx , ~, ax] = NSCA(xx, Ix, Jx);
[sy , ~, ay] = NSCA(yy, Iy, Jy);

% denoising by deflation
xxx = xx;
yyy = yy;
qtils = 0.9:-.1:.6;
for k = 1:length(qtils),
    Ex = sqrt(filter(ones(1, lx), lx, xxx.^2, [], 2));
    Ey = sqrt(filter(ones(1, ly), ly, yyy.^2, [], 2));
    
    eex = mean(Ex, 1);
    eey = mean(Ey, 1);
    
    thx = quantile(eex, qtils(k));
    thy = quantile(eey, qtils(k));
    Ix = find(eex >= thx);
    Jx = 1:length(Ex);%find(eex < thx);%
    Iy = find(eey >= thy);
    Jy = 1:length(Ey);%find(eey < thy);%
    [sxx , ~, ax] = NSCA(xxx, Ix, Jx);
    [syy , ~, ay] = NSCA(yyy, Iy, Jy);
    
    for p = 1:3,%size(sxx,1),
        sxx(p, :) = sxx(p, :) - wden(sxx(p, :), 'rigrsure', 's', 'mln', 4, 'db8');
        syy(p, :) = syy(p, :) - wden(syy(p, :), 'rigrsure', 's', 'mln', 4, 'db8');
    end
    
    xxx_old = xxx;
    yyy_old = yyy;
    
    xxx = ax*sxx;
    yyy = ay*syy;
    %     tx = (0:length(xxx)-1)/fsx;
    %     for i = 1:size(x, 1),
    %         figure;
    %         hold on;
    %         plot(tx, xxx_old(i, :), 'r');
    %         plot(tx, xxx(i, :), 'b');
    %         grid
    %         title(num2str(i));
    %     end
    %
    %     ty = (0:length(yyy)-1)/fsy;
    %     for i = 1:size(y, 1),
    %         figure;
    %         hold on;
    %         plot(ty, yyy_old(i, :), 'r');
    %         plot(ty, yyy(i, :), 'b');
    %         grid
    %         title(num2str(i));
    %     end
end
% [sx , ~, ax] = NSCA2(xx, Ix, Jx, qtl);
% [sy , ~, ay] = NSCA2(yy, Iy, Jy, qtl);

% Ix = 1:length(xx);
% Jx = 1:length(xx);
% wxI = linspace(wlow, whigh, length(Ix));
% % wxI = logspace(0, 1, length(Ix));
% wxJ = linspace(1, 1, length(Jx));
% [sx , ~, ax] = WNSCA(xx, Ix, wxI, Jx, wxJ);
% Iy = 1:length(yy);
% Jy = 1:length(yy);
% wyI = linspace(wlow, whigh, length(Iy));
% % wyI = logspace(0, 1, length(Iy));
% wyJ = linspace(1, 1, length(Jy));
% [sy , ~, ay] = WNSCA(yy, Iy, wyI, Jy, wyJ);

channelnumx = size(x, 1);
channelnumy = size(y, 1);
epochnum = floor((length(x)-epochlen*fsx)/((epochlen-epochoverlap)*fsx));
wxPCA = zeros(epochnum, sourcenum, channelnumx);
wyPCA = zeros(epochnum, sourcenum, channelnumy);
wxJADE = zeros(epochnum, sourcenum, channelnumx);
wyJADE = zeros(epochnum, sourcenum, channelnumy);
wxxJADE = zeros(epochnum, sourcenum, channelnumx);
wyyJADE = zeros(epochnum, sourcenum, channelnumy);
for i = 1:epochnum,
    xx_ = xx(:, round(i*(epochlen-epochoverlap)*fsx:i*(epochlen-epochoverlap)*fsx+epochlen*fsx));
    yy_ = yy(:, round(i*(epochlen-epochoverlap)*fsy:i*(epochlen-epochoverlap)*fsy+epochlen*fsy));
    
    % PCA
    Cx = cov(xx_');
    [V D] = eig(Cx);
    [~,I] = sort(diag(D), 1, 'descend');
    wxPCA(i,:,:) = V(:, I(1:sourcenum))';
    
    Cy = cov(yy_');
    [V D] = eig(Cy);
    [~,I] = sort(diag(D), 1, 'descend');
    wyPCA(i,:,:) = V(:, I(1:sourcenum))';

    % JADE
    wxJADE(i,:,:) = jadeR(xx_, sourcenum);
    wxxJADE(i,:,:) = normr(squeeze(wxJADE(i,:,:)));
    
    wyJADE(i,:,:) = jadeR(yy_, sourcenum);
    wyyJADE(i,:,:) = normr(squeeze(wyJADE(i,:,:)));

    % SOBI
    %             w = pinv(sobi(xx));
    %             wSOBI(i,:,:) = w(1:sourcenum, :);
    
    % SCA
    %             [~, w] = SCA2(xx, 10.0/fs, 10.0/fs, 5);
    %             wSCA(i,:,:) = w(1:sourcenum, :);
    
    %     figure
    %     plot(xx');
    %     grid
    %     title(num2str(i));
end
AnglesxPCA = zeros(epochnum, 1);
AnglesyPCA = zeros(epochnum, 1);
AnglesxJADE = zeros(epochnum, 1);
AnglesyJADE = zeros(epochnum, 1);
%         AnglesSOBI = zeros(epochnum, 1);
%         AnglesSCA = zeros(epochnum, 1);
for i = 1:epochnum, % over epochs
    AnglesxPCA(i) = subspace(squeeze(wxPCA(i,:,:))', squeeze(wxPCA(1,:,:))')*180/pi;
    AnglesyPCA(i) = subspace(squeeze(wyPCA(i,:,:))', squeeze(wyPCA(1,:,:))')*180/pi;

    AnglesxJADE(i) = subspace(squeeze(wxJADE(i,:,:))', squeeze(wxJADE(1,:,:))')*180/pi;
    AnglesyJADE(i) = subspace(squeeze(wyJADE(i,:,:))', squeeze(wyJADE(1,:,:))')*180/pi;
end
% llx = round(wlen2*fsx);
% lly = round(wlen2*fsy);
%
% Esx = sqrt(filter(ones(1, llx), llx, sx.^2, [], 2));
% Esy = sqrt(filter(ones(1, lly), lly, sy.^2, [], 2));
% Esx = zeros(size(sx));
% Esy = zeros(size(sy));
% for k = 1:size(sx, 2),
%     ind = max(k-round(llx/2), 1):min(k+round(llx/2), length(sx));
%     Esx(:, k) = median(abs(sx(:, ind)), 2);
% end
% for k = 1:size(sy, 2),
%     ind = max(k-round(lly/2), 1):min(k+round(lly/2), length(sy));
%     Esy(:, k) = median(abs(sy(:, ind)), 2);
% end

% xx = ax(:, 1:keep) * sx(1:keep, :);
% yy = ay(:, 1:keep) * sy(1:keep, :);
% xx = ax * (sx - ssx);
% yy = ay * (sy - ssy);

% IIR notch
% Wo = f0/(ffs/2);  BW = Wo/Q;
% [b,a] = iirnotch(Wo,BW);
% x = filter(b, a, x, [], 2);
% y = filter(b, a, y, [], 2);

% xxLP = LPFilter(x, 5.0/fsx);
% xxBP = BPFilter(x, 7.0/fsx, 15.0/fsx);
% yyLP = LPFilter(y, 5.0/fsy);
% yyBP = BPFilter(y, 7.0/fsy, 15.0/fsy);
%
% rx = sum(xxLP.^2, 2)./sum(xxBP.^2, 2);
% ry = sum(yyLP.^2, 2)./sum(yyBP.^2, 2);
%
% medpowerradio_x = median(rx)
% medpowerradio_y = median(ry)

% % estimate spectrum
% h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% % h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% % h = spectrum.periodogram('Rectangular');
% % h = spectrum.yulear;
% Hx = MultiChannelSpectrum(x, fsx, nfft, h); % Estimate the PSD
% Hy = MultiChannelSpectrum(y, fsy, nfft, h); % Estimate the PSD

% [x, W, A] = SCA(x, .5/fsy, 7/fsx);
% [y, W, A] = SCA(y, .5/fsy, 7/fsy);

% JADE
wx = jadeR(xx, 8);
qx = wx*xx;
wy = jadeR(yy, 8);
qy = wy*yy;

wxx = jadeR(xxx, 8);
px = wxx*xxx;
wyy = jadeR(yyy, 8);
py = wyy*yyy;
% Ex = sqrt(filter(ones(1, lx), lx, xxx.^2, [], 2));
% Ey = sqrt(filter(ones(1, ly), ly, yyy.^2, [], 2));
% eex = mean(Ex, 1);
% eey = mean(Ey, 1);
% thx = quantile(eex, probthx);
% thy = quantile(eey, probthy);
% Ix = find(eex >= thx);
% Jx = 1:length(Ex);%find(eex < thx);%
% Iy = find(eey >= thy);
% Jy = 1:length(Ey);%find(eey < thy);%
% [px , ~, ax] = NSCA(xxx, Ix, Jx);
% [py , ~, ay] = NSCA(yyy, Iy, Jy);

% SOBI
% [~, x] = sobi(x);
% [~, y] = sobi(y);


% [Vx Dx] = eig(Cx);
% [Vy Dy] = eig(Cy);
%
% Dx = diag(Dx);
% Dy = diag(Dy);
%
% Dx = Dx(end:-1:1);
% Dy = Dy(end:-1:1);
%
% Dx = Dx/Dx(1);
% Dy = Dy/Dy(1);
%
% figure
% hold on
% bar(Dx,'b');
% bar(Dy,'r');
% h = findobj(gca,'Type','patch');
% set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
% set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
% grid

% figure
% hold on
% for i = 1:size(x, 1),
%     plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
%     plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
% end
% grid;

% tx = (0:length(x)-1)/fsx;
% txx = (0:length(xx)-1)/fsx;
% for i = 1:size(x, 1),
%     figure;
%     hold on;
%     plot(tx, x(i, :), 'b');
%     plot(txx, xx(i, :), 'r');
%     grid
%     title(num2str(i));
% end
%
% ty = (0:length(y)-1)/fsy;
% tyy = (0:length(yy)-1)/fsy;
% for i = 1:size(y, 1),
%     figure;
%     hold on;
%     plot(ty, y(i, :), 'k');
%     plot(tyy, yy(i, :), 'm');
%     grid
%     title(num2str(i));
% end

% PlotECG(x, 4, 'b', fsx);
% PlotECG(y, 4, 'k', fsy);

PlotECG(qx(1:5,:), 5, 'b', fsx, 'ICA on x');
PlotECG(qy(1:5,:), 5, 'r', fsy, 'ICA on y');
%
PlotECG(sx(1:5,:), 5, 'b', fsx, 'NSCA on x');
PlotECG(sy(1:5,:), 5, 'r', fsy, 'NSCA on y');
%
PlotECG(px(1:5,:), 5, 'b', fsx, 'NSCA on x after deflation');
PlotECG(py(1:5,:), 5, 'r', fsy, 'NSCA on y after deflation');

% tx = (0:length(sx)-1)/fsx;
% for i = 1:3,%size(sx, 1),
%     figure;
%     hold on;
%     plot(tx, sx(i, :), 'b');
%     plot(tx, Esx(i, :), 'r');
%     grid
%     title(num2str(i));
% end
%
% ty = (0:length(sy)-1)/fsy;
% for i = 1:3,%size(sy, 1),
%     figure;
%     hold on;
%     plot(tx, sy(i, :), 'k');
%     plot(tx, Esy(i, :), 'm');
%     grid
%     title(num2str(i));
% end

nwindow = round(fsx*SpectralWinLen);
noverlap = round(fsx*SpectralOverlapPercentage*SpectralWinLen/100.0);
for c = 1:5,%min(size(x, 1), size(y, 1));
    figure
    subplot(211);
    [~,f,t,p] = spectrogram(sx(c, :), nwindow, noverlap, nfft, fsx, 'yaxis');
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    xlabel('Time');
    ylabel('Frequency (Hz)');
    title(num2str(c));
    subplot(212);
    [~,f,t,p] = spectrogram(sy(c, :), nwindow, noverlap, nfft, fsy, 'yaxis');
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    xlabel('Time');
    ylabel('Frequency (Hz)');
end

tx = (0:length(AnglesxJADE)-1)/fsx;
ty = (0:length(AnglesyJADE)-1)/fsy;
figure
subplot(211);
% plot(AnglesxPCA);
% hold on;
plot(tx, AnglesxJADE, 'r');
grid
title('x');
subplot(212);
% plot(AnglesyPCA);
% hold on;
plot(ty, AnglesyJADE,'r');
grid
title('y');
