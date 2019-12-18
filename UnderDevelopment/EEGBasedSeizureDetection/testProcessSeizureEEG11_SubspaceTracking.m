close all;
clear;
clc;

% parameters
subject1 = 'Dog_3';
subject2 = 'Dog_3';
number1 = 8;
number2 = 8;
path1 = ['J:\Seizure\' subject1 '\'];
path2 = ['J:\Seizure\' subject2 '\'];
mode1 = 'interictal';%'preictal';
mode2 = 'preictal';%'interictal';
epochlen = 30; % epoch length in seconds
epochoverlap = 20; % overlap between epochs in seconds
sourcenum = 3;
% nfft = 5000 - 1;
% SpectralWinLen = 10; % in seconds
% SpectralOverlapPercentage = 50; % [0 99]

% load signals
% [x fsx sequencex] = LoadSeizureEEG(path, subject, mode1, number1, ch, time);
% [y fsy sequencey] = LoadSeizureEEG(path, subject, mode2, number2, ch, time);
[x0 fsx] = LoadSeizureEEGFull(path1, subject1, mode1, number1);
[y0 fsy] = LoadSeizureEEGFull(path2, subject2, mode2, number2);

lev = 8;
x = zeros(size(x0));
for i = 1:size(x, 1)
    x(i, :) = wden(x0(i, :),'heursure','s','sln',lev,'coif5');
end

y = zeros(size(y0));
for i = 1:size(y, 1)
    y(i, :) = wden(y0(i, :),'heursure','s','sln',lev,'coif5');
end

tx = (0:size(x, 2)-1)/fsx;
for i = 1:size(x, 1),
    figure;
    hold on;
    plot(tx, x0(i, :),'b');
    plot(tx, x(i, :),'r');
    grid
    title(num2str(i));
end

ty = (0:size(y, 2)-1)/fsy;
for i = 1:size(y, 1),
    figure;
    hold on;
    plot(y0(i, :),'k');
    plot(y(i, :),'m');
    grid
    title(num2str(i));
end

% % % % % % epochnum = 230; % number of epochs
% % % % % epochnum = floor((length(x)-epochlen*fsx)/((epochlen-epochoverlap)*fsx));
% % % % % 
% % % % % % pre-process
% % % % % % % % x = LPFilter(x, 50.0/fsx); % lowpass filter
% % % % % % % % x = x - LPFilter(x, 0.1/fsx); % highpass filter
% % % % % % % % %
% % % % % % % % y = LPFilter(y, 50.0/fsy); % lowpass filter
% % % % % % % % y = y - LPFilter(y, 0.1/fsy); % highpass filter
% % % % % 
% % % % % channelnumx = size(x,1);
% % % % % channelnumy = size(y,1);
% % % % % wxx = zeros(epochnum, sourcenum, channelnumx);
% % % % % wyy = zeros(epochnum, sourcenum, channelnumy);
% % % % % Anglesx = zeros(epochnum, 1);
% % % % % Anglesy = zeros(epochnum, 1);
% % % % % Qmatrixx = zeros(epochnum, sourcenum, sourcenum);
% % % % % Qmatrixy = zeros(epochnum, sourcenum, sourcenum);
% % % % % xx1 = x(:, round((epochlen-epochoverlap)*fsx:(epochlen-epochoverlap)*fsx+epochlen*fsx));
% % % % % yy1 = y(:, round((epochlen-epochoverlap)*fsy:(epochlen-epochoverlap)*fsy+epochlen*fsy));
% % % % % for i = 1:epochnum,
% % % % %     xx = x(:, round(i*(epochlen-epochoverlap)*fsx:i*(epochlen-epochoverlap)*fsx+epochlen*fsx));
% % % % %     yy = y(:, round(i*(epochlen-epochoverlap)*fsy:i*(epochlen-epochoverlap)*fsy+epochlen*fsy));
% % % % %     
% % % % %     % JADE
% % % % %     wxx(i,:,:) = jadeR(xx, sourcenum);
% % % % %     wyy(i,:,:) = jadeR(yy, sourcenum);
% % % % %     %     wxx(i,:,:) = normr(jadeR(xx, sourcenum));
% % % % %     %     wyy(i,:,:) = normr(jadeR(yy, sourcenum));
% % % % %     
% % % % %     % SOBI
% % % % %     %         wx = pinv(sobi(xx));
% % % % %     %         wxx(i,:,:) = normr(wx(1:sourcenum, :));
% % % % %     %         wy = pinv(sobi(yy));
% % % % %     %         wyy(i,:,:) = normr(wy(1:sourcenum, :));
% % % % %     
% % % % %     % SCA
% % % % %     %     [~, wx] = SCA2(xx, 10.0/fsx, 10.0/fsx, 5);
% % % % %     %     wxx(i,:,:) = normr(wx(1:sourcenum, :));
% % % % %     %     [~, wy] = SCA2(yy, 10.0/fsy, 10.0/fsy, 5);
% % % % %     %     wyy(i,:,:) = normr(wy(1:sourcenum, :));
% % % % %     
% % % % %     %     figure
% % % % %     %     plot(xx');
% % % % %     %     grid
% % % % %     %     title(num2str(i));
% % % % %     
% % % % %     %     Anglesx(i) = subspace(xx', xx1')*180/pi;
% % % % %     %     Anglesy(i) = subspace(yy', yy1')*180/pi;
% % % % %     
% % % % %     Anglesx(i) = subspace(squeeze(wxx(i,:,:))', squeeze(wxx(1,:,:))')*180/pi;
% % % % %     Anglesy(i) = subspace(squeeze(wyy(i,:,:))', squeeze(wyy(1,:,:))')*180/pi;
% % % % %     
% % % % %     %     [Q,R] = QR(A,0);
% % % % %     %
% % % % %     %     Cx = squeeze(wxx(i,:,:))*squeeze(wxx(1,:,:))';
% % % % %     %     [Ux,Sx,Vx] = svd(Cx);
% % % % %     %     Qx = Ux*Vx';
% % % % %     %     Qx = Qx';
% % % % %     %     Qmatrixx(i,:,:) = Qx;
% % % % %     %     Anglesx(i, :) = acos(diag(squeeze(Qmatrixx(i,:,:))'*squeeze(Qmatrixx(1,:,:))))*180/pi;
% % % % %     %
% % % % %     %     Cy = squeeze(wyy(i,:,:))*squeeze(wyy(1,:,:))';
% % % % %     %     [Uy,Sy,Vy] = svd(Cy);
% % % % %     %     Qy = Uy*Vy';
% % % % %     %     Qy = Qy';
% % % % %     %     Qmatrixy(i,:,:) = Qy;
% % % % %     %     Anglesy(i, :) = acos(diag(squeeze(Qmatrixy(i,:,:))'*squeeze(Qmatrixy(1,:,:))))*180/pi;
% % % % %     
% % % % %     disp(['epoch = ', num2str(i)]);
% % % % % end
% % % % % 
% % % % % figure(5);
% % % % % subplot(211);
% % % % % plot(Anglesx, 'b');
% % % % % a = axis; a(3) = 0; a(4) = 90; axis(a);
% % % % % grid;
% % % % % title(mode1);
% % % % % subplot(212);
% % % % % plot(Anglesy, 'r');
% % % % % a = axis; a(3) = 0; a(4) = 90; axis(a);
% % % % % grid;
% % % % % title(mode2);
% % % % % 
% % % % % % pre-process
% % % % % % x = LPFilter(x, 50.0/fsx); % lowpass filter
% % % % % % x = x - LPFilter(x, 1.0/fsx); % highpass filter
% % % % % % %
% % % % % % y = LPFilter(y, 50.0/fsy); % lowpass filter
% % % % % % y = y - LPFilter(y, 1.0/fsy); % highpass filter
% % % % % 
% % % % % % xxLP = LPFilter(x, 5.0/fsx);
% % % % % % xxBP = BPFilter(x, 7.0/fsx, 15.0/fsx);
% % % % % % yyLP = LPFilter(y, 5.0/fsy);
% % % % % % yyBP = BPFilter(y, 7.0/fsy, 15.0/fsy);
% % % % % %
% % % % % % rx = sum(xxLP.^2, 2)./sum(xxBP.^2, 2);
% % % % % % ry = sum(yyLP.^2, 2)./sum(yyBP.^2, 2);
% % % % % %
% % % % % % medpowerradio_x = median(rx)
% % % % % % medpowerradio_y = median(ry)
% % % % % 
% % % % % % wx = jadeR(x);
% % % % % % wy = jadeR(y);
% % % % % % sx = wx*x;
% % % % % % sy = wy*y;
% % % % % 
% % % % % % [wxinv, sx] = sobi(x);
% % % % % % [wyinv, sy] = sobi(y);
% % % % % 
% % % % % % fl = 0.5;
% % % % % % fu = 5.0;
% % % % % % sx = SCA2(x,fl/fsx,fu/fsx);
% % % % % % sy = SCA2(y,fl/fsy,fu/fsy);
% % % % % 
% % % % % % % estimate spectrum
% % % % % % h = spectrum.welch('Hamming', round(SpectralWinLen*fsx), SpectralOverlapPercentage); % Create a Welch spectral estimator.
% % % % % % % h = spectrum.periodogram('Hamming'); % Create a Periodogram spectral estimator.
% % % % % % % h = spectrum.periodogram('Rectangular');
% % % % % % % h = spectrum.yulear;
% % % % % % Hx = MultiChannelSpectrum(sx, fsx, nfft, h); % Estimate the PSD
% % % % % % Hy = MultiChannelSpectrum(sy, fsy, nfft, h); % Estimate the PSD
% % % % % 
% % % % % % Cx = cov(x');
% % % % % % Cy = cov(y');
% % % % % %
% % % % % % [Vx Dx] = eig(Cx);
% % % % % % [Vy Dy] = eig(Cy);
% % % % % % Dx = diag(Dx);
% % % % % % Dy = diag(Dy);
% % % % % % Dx = Dx(end:-1:1);
% % % % % % Dy = Dy(end:-1:1);
% % % % % % Dx = Dx/Dx(1);
% % % % % % Dy = Dy/Dy(1);
% % % % % 
% % % % % % figure
% % % % % % hold on
% % % % % % bar(Dx,'b');
% % % % % % bar(Dy,'r');
% % % % % % h = findobj(gca,'Type','patch');
% % % % % % set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
% % % % % % set(h(1), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
% % % % % % grid
% % % % % 
% % % % % % figure
% % % % % % hold on
% % % % % % for i = 1:size(x, 1),
% % % % % %     plot(Hx(i).Frequencies, Hx(i).data/sqrt(mean(Hx(i).data.^2)),'b');
% % % % % %     plot(Hy(i).Frequencies, Hy(i).data/sqrt(mean(Hy(i).data.^2)),'r');
% % % % % % end
% % % % % % grid;
% % % % % 
% % % % % %
% % % % % % PlotECG(sx, 4, 'b', fsx);
% % % % % % PlotECG(sy, 4, 'r', fsy);
% % % % % 
% % % % % % for i = 1:size(x, 1),
% % % % % %     figure;
% % % % % %     hold on;
% % % % % %     plot(x0(i, :),'b');
% % % % % %     plot(x(i, :),'r');
% % % % % %     grid
% % % % % %     title(num2str(i));
% % % % % % end
% % % % % % 
% % % % % % for i = 1:size(y, 1),
% % % % % %     figure;
% % % % % %     hold on;
% % % % % %     plot(y0(i, :),'k');
% % % % % %     plot(y(i, :),'m');
% % % % % %     grid
% % % % % %     title(num2str(i));
% % % % % % end
% % % % % 
% % % % % % for i = 1:size(x, 1),
% % % % % %     figure;
% % % % % %     hold on;
% % % % % %     plot(cy(i, :),'r');
% % % % % %     plot(cx(i, :));
% % % % % %     grid
% % % % % %     title(num2str(i));
% % % % % % end
% % % % % 
% % % % % % c = 1;
% % % % % % nwindow = round(fsx*3.0);
% % % % % % noverlap = round(fsx*2.5);
% % % % % % nfft = 1000;
% % % % % %
% % % % % % figure
% % % % % % subplot(211);
% % % % % % spectrogram(x(c, :), nwindow, noverlap, nfft, fsx, 'yaxis');
% % % % % % subplot(212);
% % % % % % spectrogram(y(c, :), nwindow, noverlap, nfft, fsy, 'yaxis');