function [quantized, gain, bias] = approximate_quantizer(data, gain_mode, bias_mode)

switch bias_mode
    case 'zero'
        bias = 0;
    case 'mean'
        bias = mean(data(:));
    case 'median'
        bias = median(data(:));
    otherwise
        error('Undefined bias calculation mode.');
end

data = data - bias;

switch gain_mode
    case 'max'
        gain = max(abs(data(:)));
        quantized = data / gain;
    case 'max_plus_eps'
        gain = max(abs(data(:)));
        quantized = data / (gain + eps);
    otherwise
        error('Undefined gain calculation mode.');

end