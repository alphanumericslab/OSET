function [sign, exponent, mantissa] = double_to_binary(x, dtype)
% FUNCTION: double_to_binary
%
% DESCRIPTION:
% This function converts a given number from either double- or single-precision IEEE Standard 754 format to its binary representation. It extracts the sign, exponent, and mantissa components of the number and returns them as separate variables.
%
% SYNTAX:
% [sign, exponent, mantissa] = double_to_binary(x, dtype)
%
% INPUT ARGUMENTS:
% - x: The number (or vector) to be converted to binary representation.
% - dtype: Data type specification. It should be either 'double' or 'single'.
%
% OUTPUT ARGUMENTS:
% - sign: The sign component of the number (-1 for negative, +1 for positive).
% - exponent: The exponent component of the number.
% - mantissa: The mantissa component of the number.
%
% USAGE:
% To convert a double precision number 'x' to binary representation and extract its components, use:
%     [sign, exponent, mantissa] = double_to_binary(x, 'double')
%
% To convert a single precision number 'x' to binary representation and extract its components, use:
%     [sign, exponent, mantissa] = double_to_binary(x, 'single')
%
% EXAMPLES:
% % Example 1: Convert a double precision number to binary representation
% x = -3.14159;
% [sign, exponent, mantissa] = double_to_binary(x, 'double');
% disp(['Sign: ', num2str(sign)]);
% disp(['Exponent: ', num2str(exponent)]);
% disp(['Mantissa: ', num2str(mantissa)]);
%
% % Example 2: Convert a single precision number to binary representation
% y = 2.71828;
% [sign, exponent, mantissa] = double_to_binary(y, 'single');
% disp(['Sign: ', num2str(sign)]);
% disp(['Exponent: ', num2str(exponent)]);
% disp(['Mantissa: ', num2str(mantissa)]);
%
% NOTE:
% - The function assumes that the input number is a valid double or single precision floating-point number (vector).
% - It supports the IEEE Standard 754 representation for double and single precision numbers.
% - The mantissa value is normalized in the range [1, 2) for double precision and [1, 2) for single precision.
% - The exponent values are adjusted by subtracting the respective biases (1023 for double precision, 127 for single precision).
% - Any other data type specified in 'dtype' will result in an error.
%
% APPLICATIONS:
% Beyond generic CS appliications, this function can be used to analyze the
% actual number of digitization bits of a given data file from an unknown
% ADC and data conversion pipeline
%
% Reza Sameni, May 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

switch dtype
    case 'single' % IEEE Standard 754 single precision number
        % Convert the single number to binary representation
        binary_rep = dec2bin(typecast(x, 'int32'), 32);

        % Extract the sign bit, exponent, and mantissa
        sign_bit = binary_rep(:, 1);
        exponent_bits = binary_rep(:, 2:9);
        mantissa_bits = binary_rep(:, 10:end);

        bias = 127;
        demonimator = 2^23;

    case 'double' % IEEE Standard 754 double precision number
        % Convert the double number to binary representation
        binary_rep = dec2bin(typecast(x, 'int64'), 64);

        % Extract the sign bit, exponent, and mantissa
        sign_bit = binary_rep(:,1);
        exponent_bits = binary_rep(:, 2:12);
        mantissa_bits = binary_rep(:, 13:end);

        bias = 1023;
        demonimator = 2^52;

    otherwise
        error("Only 'double' and 'single' data-types are supported");

end

% Convert the bits to decimal representation
sign = (-1)^(str2double(sign_bit));
exponent = bin2dec(exponent_bits) - bias;
mantissa = bin2dec(mantissa_bits)/demonimator;
