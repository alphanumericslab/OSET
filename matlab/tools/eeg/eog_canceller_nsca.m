function [eeg_refined, Delta, delta] = eog_canceller_nsca(eeg, eog, fs, power_env_wlen, num_itr, th_method, th_value, num_comp_est_method, shrinkage_method, shrinkage_params, adapt_reference, flagplot)
%
% Implements an Electrooculography (EOG) artifact removal algorithm from
% Electroencephalography (EEG) data using Non-Stationary Component Analysis (NSCA).
% The method iteratively identifies and removes or reduces EOG-related components
% from EEG recordings based on the deflation approach.
%
% Syntax:
% [eeg_refined, Delta, delta] = eog_canceller_nsca(eeg, eog, fs, power_env_wlen, num_itr, th_method, th_value, num_comp_est_method, shrinkage_method, shrinkage_params, adapt_reference, flagplot)
%
% Inputs:
%   eeg: a [channels x samples] matrix of noisy EEG data.
%   eog: a [channels x samples] matrix or vector of EOG signals. Specify as [] to use the first EEG channel as EOG (not recommended).
%   fs: the sampling frequency in Hz.
%   power_env_wlen: the half-length of the window (in seconds) for calculating the power envelope.
%   num_itr: maximum number of iterations for the artifact removal process.
%   th_method: method for selecting the threshold for EOG segment detection; options include 'absolute', 'envelope-mean-fraction', 'envelope-median-fraction', 'envelope-prctile'.
%   th_value: numeric value specifying the threshold level, whose interpretation depends on th_method.
%   num_comp_est_method: method to estimate the number of components to denoise each iteration ('fixed' number defined by shrinkage_params.num_noise_gevds, or 'opt' for optimal estimation).
%   shrinkage_method: technique for shrinking EOG components. 'remove' for setting to zero (reduces dimension), 'wden' (for wavelet denoising).
%   shrinkage_params: structure or parameters specific to the shrinkage method selected.
%   adapt_reference: controls if the reference EOG should be adapted across iterations ('no', 'top-gevd', 'top-few-gevd').
%   flagplot: specifies the verbosity of the plots (0 = none, 1 or 2 for increasing detail).
%
% Outputs:
%   eeg_refined: The filtered EEG data as a [channels x samples] matrix.
%   Delta: A vector containing the relative power of EEG data in high-power segments for each iteration.
%   delta: A matrix with rows corresponding to iterations and columns to channels, showing relative power per channel in high-power segments.
%
% Notes:
%   - It is highly recommended to use actual EOG channels for better artifact removal instead of leaving the EOG field empty.
%   - The function's performance and effectiveness depend on sampling frequency and appropriate parameter selection, which might require experimentation based on the specific dataset.
%   - See dimension_estimation() and reference below for num_comp_est_method = 'opt' (optimal estimation)
%
% References:
%   Sameni, R., and Gouy-Pailler, C. (2014). An iterative subspace denoising
%   algorithm for removing electroencephalogram ocular artifacts.
%   In Journal of Neuroscience Methods (Vol. 225, pp. 97–105).
%   https://doi.org/10.1016/j.jneumeth.2014.01.024
% 
%   Sameni, R., Jutten, C., & Shamsollahi, M. B. (2010). A Deflation
%   Procedure for Subspace Decomposition. In IEEE Transactions on Signal
%   Processing (Vol. 58, Issue 4, pp. 2363–2374). Institute of Electrical
%   and Electronics Engineers (IEEE).
%   https://doi.org/10.1109/tsp.2009.2037353
%
% Revision History:
%   2008: First release.
%   2024: Renamed from deprecated version EOGRemoval().
%       Comprehensive documentation update for improved clarity and usability.
%
% Developed by Reza Sameni, 2008-2024
% The Open-Source Electrophysiological Toolbox (OSET)

%

if isequal(shrinkage_method, 'remove') && num_itr > 1
    error('when shrinkage_method = ''remove'', num_itr must be 1 (single step)')
end

if isempty(eog)
    eog = eeg(1,:);
end

% define the parameters
power_env_wlen_in_samples = round(fs * power_env_wlen);
num_ch = size(eeg, 1);
num_samp = size(eeg, 2);
reference_channel = sqrt(mean(eog.^2, 1));
eeg_refined = eeg;

Delta = [];
delta = [];
for i = 1 : num_itr
    % calculate the time-varying power envelope
    power_envelope = sqrt(filtfilt(ones(1,power_env_wlen_in_samples), power_env_wlen_in_samples, reference_channel.^2)); % notice the normalization by SD of EOG, and not the SD of reference_channel!

    % threshold the power envelope
    switch th_method
        case 'absolute'
            I_high_power = power_envelope > th_value;
            I_low_power = power_envelope <= th_value;
        case 'envelope-mean-fraction'
            I_high_power = power_envelope > th_value * mean(power_envelope);
            I_low_power = power_envelope <= th_value * mean(power_envelope);
        case 'envelope-median-fraction'
            I_high_power = power_envelope > th_value * median(power_envelope);
            I_low_power = power_envelope <= th_value * median(power_envelope);
        case 'envelope-prctile'
            th = prctile(power_envelope, th_value);
            I_high_power = power_envelope > th;
            I_low_power = power_envelope <= th;
        otherwise
            error('undefined power thresholding method');
    end
    % plot the detected EOG segments
    if flagplot == 1
        p = zeros(1,num_samp);
        p(I_high_power) = 1;
        t = (0:num_samp-1)/fs;
        figure;
        subplot(311);
        plot(t, reference_channel,'k','linewidth',0.5);
        grid;
        set(gca,'Box','On','FontSize',16);
        set(gca,'XTickLabel',[]);
        ylabel('Amplitude(mv)');
        axis tight;

        subplot(312);
        plot(t, power_envelope,'k','linewidth',2);
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
    elseif flagplot == 2
        t = (0:num_samp-1)/fs;
        figure;
        plot(t, reference_channel,'b','linewidth',0.5);
        hold on
        plot(t, power_envelope,'r','linewidth',2);
        plot(t(I_high_power),power_envelope(I_high_power),'k.','linewidth',2);
        grid;
        set(gca,'Box','On','FontSize',16);
        legend('Reference channel', 'Power envelope', 'Detected segments');
        ylabel('Amplitude(mv)');
        xlabel('time(s)');
        axis tight;
    end

    if isempty(find(I_high_power,1)) || isempty(find(I_low_power,1))
        break
    end

    % calculate the delta and Delta parameters
    Delta0 = trace(eeg_refined(:,I_high_power)*eeg_refined(:,I_high_power)')/trace(eeg_refined(:,1:num_samp)*eeg_refined(:,1:num_samp)');
    delta0 = diag(eeg_refined(:,I_high_power)*eeg_refined(:,I_high_power)')./diag(eeg_refined(:,1:num_samp)*eeg_refined(:,1:num_samp)');

    Delta = cat(1, Delta, Delta0);
    delta = cat(1, delta, delta0);

    % nonstationary component analysis (GEVD over given windows)
    [y, ~, A] = nonstationary_component_analysis(eeg_refined, I_high_power, 1:num_samp);
    A = real(A);
    y = real(y);

    switch num_comp_est_method
        case 'fixed' % fixed number of channels corresponding to the largest generalized eigenvalues
            num_ch_den_per_itr = shrinkage_params.num_noise_gevds;
        case 'opt' % estimate the EOG dimensions embedded in background EEG noise
            eeg_refined_normalized = (eeg_refined - mean(eeg_refined,2))./repmat(std(eeg_refined(:,I_low_power),[],2), 1, num_samp); % normalize data
            [~, ~, ~, ~, ~, num_ch_den_per_itr] = dimension_estimation(eeg_refined_normalized(:,I_high_power), var(eeg_refined_normalized(1,I_low_power)), flagplot == 3);
        otherwise
            error('undefined shrinkage method');
    end
    if num_ch_den_per_itr >= num_ch
        warning('num_ch_den_per_itr more than or equal to the number of channels');
        num_ch_den_per_itr = num_ch - 1;
    end

    switch shrinkage_method
        case 'remove'
            eeg_refined = A(:, num_ch_den_per_itr + 1 : end) * y(num_ch_den_per_itr + 1 : end, :);
        case 'wden' % wavelet denoising
            y_den = y;
            for k = 1 : num_ch_den_per_itr
                eog_est = wden(y(k,:), shrinkage_params.TPTR, shrinkage_params.SORH, shrinkage_params.SCAL, shrinkage_params.NUM, shrinkage_params.WNAME);
                y_den(k,:) = y_den(k,:) - eog_est;
            end
            eeg_refined = A * y_den;
        otherwise
            error('undefined shrinkage method');
    end

    % recalculate the reference channel
    switch adapt_reference
        case 'no'
            % keep the reference
        case 'top-gevd'
            reference_channel = y(1,:);
        case 'top-few-gevd'
            reference_channel = sqrt(mean(y(1:num_ch_den_per_itr,:).^2, 1));
        otherwise
            error('Undefined reference adaptation method');
    end

end