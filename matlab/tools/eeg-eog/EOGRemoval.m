function x = EOGRemoval(EEG,EOG,wlen,th,M,L,adapt,flagplot)
%
% x = EOGRemoval(EEG,EOG,wlen,th,M,L,adapt,flagplot),
% EOG removal from EEG by deflation
%
% inputs:
% EEG: matrix of noisy data (channels x samples)
% EOG: vector of ocular signal (1 x samples). By default EEG(1,:) is used
% if EOG is an empty vector; but don't do that please! :)
% wlen: the half length of the energy window (in samples)
% th: the energy threshold for finding active EOG times
% M: number of channels to denoise in each iteration
% L: number of iterations (the algorithm may stop before reaching L)
% adapt: to adapt(1) or leave unchanged(0) the reference EOG
% flagplot: to plot(1) or not to plot(0) the results
%
% output:
% b: vector or matrix of baseline wanders (channels x samples)
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

if(isempty(EOG))
    EOG = EEG(1,:);
end

% define the parameters


N = size(EEG,2);
ref = EOG;
x = EEG;
% energy = zeros(1,N);

typical = x(1,:);
Delta = [];
delta = [];
for i = 1:L
    % calculate the time-varying energy
    energy = sqrt(filtfilt(ones(1,wlen),wlen,ref.^2))./std(EOG); % notice the normalization by SD of EOG, and not the SD of ref!

    % threshold the energy
    I = energy > th*mean(energy);
    J = energy <= th*mean(energy);
    if(isempty(find(I,1)))
        break;
    end

    % calculate the delta and Delta parameters
    Delta0 = trace(x(:,I)*x(:,I)')/trace(x(:,1:N)*x(:,1:N)');
    Delta0 = round(1000*Delta0)/10;
    delta0 = diag(x(:,I)*x(:,I)')./diag(x(:,1:N)*x(:,1:N)');
    delta0 = round(1000*delta0)/10;

    p = zeros(1,N);
    p(I) = 1;

    fs = 250;
    t = (0:N-1)/fs;

    if(flagplot==2)
        %     plot the results
        figure;
        subplot(311);
        plot(t,ref,'k','linewidth',0.5);
        grid;
        set(gca,'Box','On','FontSize',16);
        set(gca,'XTickLabel',[]);
        ylabel('Amplitude(mv)');
        axis tight;

        subplot(312);
        plot(t,energy,'k','linewidth',2);
        grid;
        set(gca,'Box','On','FontSize',16);
        set(gca,'XTickLabel',[]);
        ylabel('Normalized Energy');

        subplot(313);
        plot(t,p,'k','linewidth',2);
        grid;
        set(gca,'Box','On','FontSize',16);
        ylabel('Activation Pulse');
        xlabel('time(s)');
    end

    % nonstationary component analysis (GEVD over given windows)
    [y, W, A, B, C] = nsca_source_separation(x,I,1:N);

    % estimate the EOG dimensions embedded in background EEG noise
    xx = (x - mean(x,2)*ones(1,size(x,2)))./(std(x(:,J),[],2)*ones(1,size(x,2))); % normalize data
    [lambda, AIC, MDL, NEW, ENSTh, ENS] = dimension_estimation(xx(:,I), var(xx(1,J)), 1);
%     ENS

    % recalculate the reference channel
    if(adapt == 1)
        ref = sqrt(mean(y(1:M,:).^2,1));
    elseif(adapt == 2)
        ref = y(1,:);
    end

    % wavelet denoising
    % % %     for k = 1:M,
    for k = 1:ENS
        if(flagplot)
            figure;
            plot(y(k,:));
            hold on;
        end
        est = wden(y(k,:),'heursure','s','mln',5,'sym5');
        % % %         est = wden(y(k,:),'rigrsure','s','mln',5,'db5');
        y(k,:) = y(k,:) - est;

        if(flagplot)
            plot(est,'g');
            grid
        end
    end

    x = A*y;

    Delta = [Delta Delta0];
    delta = [delta delta0];

    typical = [typical ; x(1,:)];
end

if(flagplot)
    % % % PlotECG(typical(:,end-fs*10+1:end),4,'b',fs);
    % % % PlotECG(typical(:,:),4,'b',fs);
    L1 = size(typical,1);
    figure;
    for i = 1:L1
        subplot(L1,1,i);
        plot(t,typical(i,:),'k','linewidth',1);
        axis([0 45 -100 100]);
        grid;
        set(gca,'Box','On','FontSize',16);
        if (i<L1)
            set(gca,'XTickLabel',[]);
        end
        %     ylabel(['IC_',num2str(i)],'FontSize',16);
        % % %     ylabel('Amplitude(mv)');
        if (i==1)
            ylabel('raw');
        else
            ylabel(num2str(i-1));
        end
    end
    xlabel('time(s)','FontSize',16);
end