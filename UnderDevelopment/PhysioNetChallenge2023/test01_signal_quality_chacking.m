clear
clc
close all

data_struct = load('../../../../DataFiles/PhysioNetChallenges2023/Sample_files_new/BIDMC_updated/ICARE_0097_20000102_185034.mat');
% data_struct = load('../../../../DataFiles/PhysioNetChallenges2023/Sample_files_new/ULB_updated/ICARE_0066_20000101_211200.mat');

fs = data_struct.sampling_rate;
data = data_struct.signal;
f_mains = 60.0;


% data = data(:, 1 : 250000);

quant_bits = 16;
max_quant = 2^(quant_bits-1) - 1;
data_quantized = int16(zeros(size(data)));
ch_gain = zeros(1, size(data, 1));
for ch = 1 : size(data, 1)
    ch_gain(ch) = max(abs(data(ch, :)));
    data_quantized(ch, :) = int16(round(max_quant * data(ch, :) / ch_gain(ch)));
end

ch_error = zeros(1, size(data, 1));
ch_snr = zeros(1, size(data, 1));
for ch = 1 : size(data, 1)
    ch_error = data(ch, :) - ch_gain(ch)*double(data_quantized(ch, :))/max_quant;
    ch_snr(ch) = 10*log10(mean(data(ch, :).^2) / mean(ch_error.^2));
end
ch_snr

num_quantization_levels = zeros(1, size(data, 1));
for ch = 1 : size(data, 1)
    num_quantization_levels(ch) = length(find(diff(sort(data(ch, :))))) + 1;
end

num_quantization_levels


x = data(:);


[sign, exponent, mantissa] = double_to_binary(x, 'double');
max_exp = max(exponent)
min_exp = min(exponent)

% [mantissa, exponent] = mantissaexp(x);
% Extract mantissa and exponent
% exponent = floor(log2(abs(x)));
% mantissa = x ./ (2.^exponent);


% I_pos = mantissa > 0;
% I_neg = mantissa < 0;
% mantissa(I_pos) = mantissa(I_pos) - 1;
% mantissa(I_neg) = mantissa(I_neg) + 1;


exponent_bit_len = ceil(log2(length(find(diff(sort(exponent)))) + 1))
mantissa_bit_len = ceil(log2(length(find(diff(sort(mantissa)))) + 1)) % - 1

figure
plot(mantissa, exponent, 'b.')
grid
xlabel('mantissa');
ylabel('exponent');
set(gca, 'fontsize', 18)


data_saturated = zeros(size(data));
k_sigma = 10;
for kk = 1 : size(data, 1)
    gamma = k_sigma * std(data(kk, :));
    data_saturated(kk, :) = gamma * tanh(data(kk, :) / gamma);
end

% Removing the mains frequency
Q = 35; % Notch filter Q-factor
w0 = f_mains/(fs/2);
BW = w0/35;
[b,a] = iirnotch(w0, BW);
data_mains_removed = filter(b, a, data_saturated')';


PlotECG(data_mains_removed, 4, 'b', fs);

t = (0:size(data, 2)-1)/fs;
for ch = 1 : size(data, 1)
    lgnd = {};
    figure;
    plot(t, data(ch, :)); lgnd = cat(1, lgnd, 'data');
    hold on
    plot(t, data_saturated(ch, :)); lgnd = cat(1, lgnd, 'data_saturated');
%     plot(t, data_mains_removed(ch, :)); lgnd = cat(1, lgnd, 'data_mains_removed');
    grid
    title(['Ch = ', num2str(ch)])
    legend(lgnd, 'interpreter', 'none')
    set(gca, 'fontsize', 18)
end


% x = data(:);
%
% y = length(find(diff(sort(x)))) + 1

