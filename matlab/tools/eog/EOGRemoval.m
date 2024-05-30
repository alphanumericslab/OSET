function x = EOGRemoval(EEG, EOG, wlen, th, M, L, adapt_reference, flagplot)
%
% x = EOGRemoval(EEG,EOG,wlen,th,M,L,adapt_reference,flagplot),
% EOG removal from EEG by deflation
%
% inputs:
% EEG: matrix of noisy data (channels x samples)
% EOG: vector of ocular signal (1 x samples). By default EEG(1,:) is used
% if EOG is an empty vector; but don't do that please! :)
% wlen: the half length of the power_envelope window (in samples)
% th: the power_envelope threshold for finding active EOG times
% M: number of channels to denoise in each iteration
% L: number of iterations (the algorithm may stop before reaching L)
% adapt_reference: to adapt_reference(1) or leave unchanged(0) the reference EOG
% flagplot: to plot(1) or not to plot(0) the results
%
% output:
% b: vector or matrix of baseline wanders (channels x samples)
%
%
%   Sameni, R. and Gouy-Pailler, C. (2014). An iterative subspace denoising
%       algorithm for removing electroencephalogram ocular artifacts.
%       In Journal of Neuroscience Methods (Vol. 225, pp. 97–105). 
%       https://doi.org/10.1016/j.jneumeth.2014.01.024
% 
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version EOGRemoval()
% 
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if isempty(EOG)
    EOG = EEG(1,:);
end

% define the parameters


N = size(EEG,2);
reference_channel = EOG;
x = EEG;
% power_envelope = zeros(1,N);

typical = x(1,:);
Delta = [];
delta = [];
for i = 1:L
    % calculate the time-varying power envelope
    power_envelope = sqrt(filtfilt(ones(1,wlen),wlen,reference_channel.^2))./std(EOG); % notice the normalization by SD of EOG, and not the SD of reference_channel!

    % threshold the power envelope
    I_high_power = power_envelope > th*mean(power_envelope);
    I_low_power = power_envelope <= th*mean(power_envelope);
    if isempty(find(I_high_power,1))
        break
    end

    % calculate the delta and Delta parameters
    Delta0 = trace(x(:,I_high_power)*x(:,I_high_power)')/trace(x(:,1:N)*x(:,1:N)');
    Delta0 = round(1000*Delta0)/10;
    delta0 = diag(x(:,I_high_power)*x(:,I_high_power)')./diag(x(:,1:N)*x(:,1:N)');
    delta0 = round(1000*delta0)/10;

    p = zeros(1,N);
    p(I_high_power) = 1;

    fs = 250;
    t = (0:N-1)/fs;

    if flagplot == 2
        %     plot the results
        figure;
        subplot(311);
        plot(t,reference_channel,'k','linewidth',0.5);
        grid;
        set(gca,'Box','On','FontSize',16);
        set(gca,'XTickLabel',[]);
        ylabel('Amplitude(mv)');
        axis tight;

        subplot(312);
        plot(t,power_envelope,'k','linewidth',2);
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
    [y, W, A, B, C] = nsca_source_separation(x, I_high_power, 1:N);

    % estimate the EOG dimensions embedded in background EEG noise
    xx = (x - mean(x,2)*ones(1,size(x,2)))./(std(x(:,I_low_power),[],2)*ones(1,size(x,2))); % normalize data
    [lambda, AIC, MDL, NEW, ENSTh, ENS] = dimension_estimation(xx(:,I_high_power), var(xx(1,I_low_power)), 1);
    %     ENS

    % recalculate the reference channel
    switch adapt_reference
        case 0
            % keep the reference
        case 1
            reference_channel = sqrt(mean(y(1:M,:).^2,1));
        case 2
            reference_channel = y(1,:);
        otherwise
            error('Undefined reference adaptation method');
    end

    % wavelet denoising
    % % %     for k = 1:M,
    for k = 1:ENS
        if flagplot
            figure;
            plot(y(k,:));
            hold on;
        end
        est = wden(y(k,:),'heursure','s','mln',5,'sym5');
        % % %         est = wden(y(k,:),'rigrsure','s','mln',5,'db5');
        y(k,:) = y(k,:) - est;

        if flagplot
            plot(est,'g');
            grid
        end
    end

    x = A*y;

    Delta = [Delta Delta0];
    delta = [delta delta0];

    typical = [typical ; x(1,:)];
end

if flagplot
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