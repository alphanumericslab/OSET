function next_state = markov_model_next_state_gen(STM, current_state)
% markov_model_next_state_gen - Generate next state from current state using a state-transition matrix for first-order Markov model.
%
% This function generates the next state from the current state based on a
% provided state-transition matrix (STM). Repeated calls to this function
% will result in state sequences with state transition probabilities
% defined by the STM.
%
% Usage Example:
%   % Define a state-transition matrix (STM)
%   STM = [0 1; 1 0];
%
%   % Initialize the first state
%   state = 1;
%
%   % Generate and display a sequence of next states
%   while 1
%       state = markov_model_next_state_gen(STM, state);
%       disp(state);
%   end
%
% Inputs:
%   STM: The state transition matrix, where STM(i,j) represents the
%       probability of transitioning from state i to state j in the next beat.
%       Each row of STM should sum up to 1.
%   current_state: The current state from which the next state is generated.
%
% Outputs:
%   next_state: The generated next state based on the STM.
%
% Notes:
%   - For each entry of STM, STM(i,j) represents the probability of
%     transitioning from state i to state j in the next beat.
%   - Each row of STM should sum up to 1.
%   - STM can be asymmetric, allowing different probabilities for
%     transitioning between states.
%
% References:
%   - Clifford, G. D., Nemati, S., & Sameni, R. (2010). An artificial vector
%     model for generating abnormal electrocardiographic rhythms. Physiological
%     Measurement, 31(5), 595–609.
%
% Revision History:
%   2023: First release

% Calculate the cumulative sum of the state transition matrix
CSTM = cumsum(STM, 2);

% Estimate the new state using random transition
a = rand;

% Determine the next state based on the cumulative probabilities
if a < CSTM(current_state, 1)
    next_state = 1;
else
    for k = 2 : size(STM, 1) % Number of states
        if a >= CSTM(current_state, k - 1) && a < CSTM(current_state, k)
            next_state = k;
            break;
        end
    end
end
