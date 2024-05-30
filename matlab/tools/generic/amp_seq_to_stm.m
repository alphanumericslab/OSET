function stm = amp_seq_to_stm(amplitude_sequence)
% amp_seq_to_stm - Generates a state transition matrix (STM) from amplitude 
%   fluctuations of a time-series. Useful for detecting periodic QRS amplitude
%   fluctuations or T-wave Alternans (TWA) from ECG data.
% 
% Usage:
%   stm = amp_seq_to_stm(amplitude_sequence)
% 
% Inputs:
%   amplitude_sequence - A numeric vector representing the amplitude sequence
%                        from a time-series data set, typically ECG data.
%
% Outputs:
%   stm - 2x2 state transition matrix that represents the transitions
%         between high and low states in the amplitude sequence.
%
%       stm(1,1): Low-to-Low transitions
%       stm(1,2): Low-to-High transitions
%       stm(2,1): High-to-Low transitions
%       stm(2,2): High-to-High transitions
% Note:
%   output STM can be converted to an estimate of high-low amplitude 
%   probabilities by the following normalization:
%   STM = diag(1./sum(stm, 2))*stm;
% 
% Example 1: White noise
%   stm = amp_seq_to_stm(randn(1, 10000)); % STM converges to [1/3, 2/3 ; 2/3, 1/3] for long
% sequences
% 
% Example 2: Wiener process
%   stm = amp_seq_to_stm(cumsum(randn(1, 10000))); % STM converges to [0.5, 0.5 ; 0.5, 0.5] for long
% sequences
% 
% Example 3: Periodic process
%   stm = amp_seq_to_stm(mod(1:10000, 2)); % STM is [0, 1 ; 1, 0]
% 
% Reference use case for T-wave alternans (TWA):
%   - Clifford, G. D., Nemati, S., & Sameni, R. (2010). An artificial vector
%     model for generating abnormal electrocardiographic rhythms. Physiological
%     Measurement, 31(5), 595–609.
%
%   Revision History:
%       2023: First release
%
%   Reza Sameni, 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Validate the input amplitude_sequence
if ~isvector(amplitude_sequence) || ~isnumeric(amplitude_sequence)
    error('input amplitude sequence must be a vector of numbers')
end

if length(amplitude_sequence) < 3
    error('input amplitude sequence must have at least 3 elements');
end

% Convert amplitude values to a binary high-low sequence.
% A value is 'high' (1) if it is greater than or equal to the preceding value, 
% and 'low' (0) otherwise.
fluctuation_sequence = amplitude_sequence(2 : end) >= amplitude_sequence(1 : end-1);

% Initialize a 2x2 state transition matrix to store the counts of transitions.
stm = zeros(2);

% Compute the number of transitions for each state.
stm(1, 1) = sum(fluctuation_sequence(1:end-1) == 0 & fluctuation_sequence(2:end) == 0); % Low to Low
stm(1, 2) = sum(fluctuation_sequence(1:end-1) == 0 & fluctuation_sequence(2:end) == 1); % Low to High
stm(2, 1) = sum(fluctuation_sequence(1:end-1) == 1 & fluctuation_sequence(2:end) == 0); % High to Low
stm(2, 2) = sum(fluctuation_sequence(1:end-1) == 1 & fluctuation_sequence(2:end) == 1); % High to High

end
