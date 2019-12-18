close all;
clear;
clc;

% parameters
subject1 = 'Dog_3';
subject2 = 'Dog_3';
number1 = 6;
number2 = 6;
path1 = ['J:\Seizure\' subject1 '\'];
path2 = ['J:\Seizure\' subject2 '\'];
mode1 = 'interictal';%'preictal';
mode2 = 'preictal';%'interictal';
epochlen = 10; % epoch length in seconds
epochoverlap = 7; % overlap between epochs in seconds
sourcenum = 3;
% nfft = 5000 - 1;
% SpectralWinLen = 10; % in seconds
% SpectralOverlapPercentage = 50; % [0 99]

% load signals
% [x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number1, ch, time);
% [y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number2, ch, time);
[x fsx] = LoadSeizureEEGFull(path1, subject1, mode1, number1);
[y fsy] = LoadSeizureEEGFull(path2, subject2, mode2, number2);

% epochnum = 230; % number of epochs
epochnum = floor((length(x)-epochlen*fsx)/((epochlen-epochoverlap)*fsx));

% pre-process
x = LPFilter(x, 10.0/fsx); % lowpass filter
x = x - LPFilter(x, 0.1/fsx); % highpass filter

y = LPFilter(y, 10.0/fsy); % lowpass filter
y = y - LPFilter(y, 0.1/fsy); % highpass filter

channelnum = size(x,1);
wxx = zeros(epochnum, sourcenum, channelnum);
wyy = zeros(epochnum, sourcenum, channelnum);
for i = 1:epochnum,
    xx = x(:, round(i*(epochlen-epochoverlap)*fsx:i*(epochlen-epochoverlap)*fsx+epochlen*fsx));
    yy = y(:, round(i*(epochlen-epochoverlap)*fsy:i*(epochlen-epochoverlap)*fsy+epochlen*fsy));
    
    % JADE
    %     wxx(i,:,:) = jadeR(xx, sourcenum);
    %     wyy(i,:,:) = jadeR(yy, sourcenum);
    wxx(i,:,:) = normr(jadeR(xx, sourcenum));
    wyy(i,:,:) = normr(jadeR(yy, sourcenum));
    
    % SOBI
    %         wx = pinv(sobi(xx));
    %         wxx(i,:,:) = normr(wx(1:sourcenum, :));
    %         wy = pinv(sobi(yy));
    %         wyy(i,:,:) = normr(wy(1:sourcenum, :));
    
    % SCA
    %     [~, wx] = SCA2(xx, 10.0/fsx, 10.0/fsx, 5);
    %     wxx(i,:,:) = normr(wx(1:sourcenum, :));
    %     [~, wy] = SCA2(yy, 10.0/fsy, 10.0/fsy, 5);
    %     wyy(i,:,:) = normr(wy(1:sourcenum, :));
    
    %     figure
    %     plot(xx');
    %     grid
    %     title(num2str(i));
    disp(['epoch = ', num2str(i)]);
end

prodsx = zeros(epochnum, epochnum, sourcenum, sourcenum);
prodsy = zeros(epochnum, epochnum, sourcenum, sourcenum);

prods2x2x = zeros(epochnum*sourcenum, epochnum*sourcenum);
prods2x2y = zeros(epochnum*sourcenum, epochnum*sourcenum);

averageonepochsx = zeros(epochnum*sourcenum, 1);
averageonepochsy = zeros(epochnum*sourcenum, 1);

epochcoherencex = zeros(epochnum, epochnum);
epochcoherencey = zeros(epochnum, epochnum);

for i = 1:epochnum, % over epochs
    for j = 1:epochnum, % over epochs
        prodsx(i, j, :, :) = abs(squeeze(wxx(i,:,:)) * squeeze(wxx(j,:,:))');
        prodsy(i, j, :, :) = abs(squeeze(wyy(i,:,:)) * squeeze(wyy(j,:,:))');
        prx = squeeze(prodsx(i, j, :, :));
        pry = squeeze(prodsy(i, j, :, :));
        
        epochcoherencex(i, j) = mean(mean(prx));
        epochcoherencey(i, j) = mean(mean(pry));
        
        %         sm = sum(pr, 2);
        %         [~, I] = sort(sm, 1, 'descend');
        %         pr = pr(I, I);
        
        prods2x2x((i-1)*sourcenum + 1 : i*sourcenum, (j-1)*sourcenum + 1 : j*sourcenum) = prx;
        prods2x2y((i-1)*sourcenum + 1 : i*sourcenum, (j-1)*sourcenum + 1 : j*sourcenum) = pry;
        
        if(i~=j)
            % % %              averageonepochs((i-1)*sourcenum + 1 : i*sourcenum, :) = averageonepochs((i-1)*sourcenum + 1 : i*sourcenum, :) + pr;
            colmaxx = max(prx, [], 2);
            colmaxy = max(pry, [], 2);
            averageonepochsx((i-1)*sourcenum + 1 : i*sourcenum, :) = averageonepochsx((i-1)*sourcenum + 1 : i*sourcenum, :) + colmaxx;
            averageonepochsy((i-1)*sourcenum + 1 : i*sourcenum, :) = averageonepochsy((i-1)*sourcenum + 1 : i*sourcenum, :) + colmaxy;
        end
    end
end

averageonepochsx = averageonepochsx/(epochnum - 1);
averageonepochsy = averageonepochsy/(epochnum - 1);

figure(1)
subplot(121);
imshow(epochcoherencex);
title('interictal');
subplot(122);
imshow(epochcoherencey);
title('preictal');

disp(['median epochcoherencex = ' num2str(median(median(epochcoherencex)))]);
disp(['median epochcoherencey = ' num2str(median(median(epochcoherencey)))]);

disp(['median prods2x2x = ' num2str(median(prods2x2x(:)))]);
disp(['median prods2x2y = ' num2str(median(prods2x2y(:)))]);

disp(['mean averageonepochsx = ' num2str(median(averageonepochsx(:)))]);
disp(['mean averageonepochsy = ' num2str(median(averageonepochsy(:)))]);

figure(2)
subplot(121);
imshow(prods2x2x);
title('interictal');
subplot(122);
imshow(prods2x2y);
title('preictal');

figure(3)
hold on
plot(averageonepochsx);
plot(averageonepochsy,'r');
grid
legend('interictal', 'preictal');

% figure(4)
% hold on
% plot(diag(epochcoherencex));
% plot(diag(epochcoherencey),'r');
% grid
% legend('interictal', 'preictal');

% PlotECG(x, 4, 'k', fsx);
% PlotECG(y, 4, 'm', fsy);

% subplot(121);
% bar(averageonepochsx);
% title('interictal');
% grid
% subplot(122);
% bar(averageonepochsy);
% title('preictal');
% grid

% pre-process
% x = LPFilter(x, 50.0/fsx); % lowpass filter
% x = x - LPFilter(x, 1.0/fsx); % highpass filter
% %
% y = LPFilter(y, 50.0/fsy); % lowpass filter
% y = y - LPFilter(y, 1.0/fsy); % highpass filter

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

% wx = jadeR(x);
% wy = jadeR(y);
% sx = wx*x;
% sy = wy*y;

% [wxinv, sx] = sobi(x);
% [wyinv, sy] = sobi(y);

% fl = 0.5;
% fu = 5.0;
% sx = SCA2(x,fl/fsx,fu/fsx);
% sy = SCA2(y,fl/fsy,fu/fsy);

% % estimate spectrum
% h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% % h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% % h = spectrum.periodogram('Rectangular');
% % h = spectrum.yulear;
% Hx = MultiChannelSpectrum(sx, fsx, nfft, h); % Estimate the PSD
% Hy = MultiChannelSpectrum(sy, fsy, nfft, h); % Estimate the PSD

% Cx = cov(x');
% Cy = cov(y');
%
% [Vx Dx] = eig(Cx);
% [Vy Dy] = eig(Cy);
% Dx = diag(Dx);
% Dy = diag(Dy);
% Dx = Dx(end:-1:1);
% Dy = Dy(end:-1:1);
% Dx = Dx/Dx(1);
% Dy = Dy/Dy(1);

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

%
% PlotECG(sx, 4, 'b', fsx);
% PlotECG(sy, 4, 'r', fsy);

% for i = 1:size(x, 1),
%     figure;
%     hold on;
%     plot(x(i, :));
%     plot(y(i, :),'r');
%     grid
%     title(num2str(i));
% end

% for i = 1:size(x, 1),
%     figure;
%     hold on;
%     plot(cy(i, :),'r');
%     plot(cx(i, :));
%     grid
%     title(num2str(i));
% end

% c = 1;
% nwindow = round(fsx*3.0);
% noverlap = round(fsx*2.5);
% nfft = 1000;
%
% figure
% subplot(211);
% spectrogram(x(c, :), nwindow, noverlap, nfft, fsx, 'yaxis');
% subplot(212);
% spectrogram(y(c, :), nwindow, noverlap, nfft, fsy, 'yaxis');