
function [P_observ_cond_to_state,P_observ_cond_to_state_comp] = gmm_pdf( P_all , mu_all  , sigma_all  , observations_time_series   ,  state_numbers ,  dimension_numbers_index , num_gmm_component  )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function work for diogonal covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_Observ_cond_to_State( t , State ) is P(y_t | x_t= State) & P_Observ_cond_to_State is a matrix (K x T)
% z = mvnpdf(X,MU,SIGMA)

%GMM parameter include GMMParam(state).P , GMMParam(state).mu , GMMParam(state).Sigma

% StateNumber is between 1 to cardinality of hidden states

% global P_Observ_cond_to_State GMMParam StateNumbers Observations nan_index correct_index miss_index P_Observ_cond_to_State_Comp

indexvalid =[];

C = length(state_numbers);
[~ , T_all]=   size(observations_time_series) ;

P_observ_cond_to_state_comp = zeros(sum(state_numbers), max(num_gmm_component), T_all);

for zee =  1:C
    current_observations = observations_time_series( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1), :);
    nan_indexes = isnan(current_observations');
    nans_patterns = unique(nan_indexes,'rows');
    
    for s = 1:state_numbers(zee)
        for k = 1:num_gmm_component(zee)
            current_mu = squeeze(mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  s , k));
            current_sig =sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1), 1:(dimension_numbers_index(zee+1)-dimension_numbers_index(zee)), s , k);
            
            for i = 1 : size (nans_patterns ,1)
                
                temp_pattern = nans_patterns(i ,:) ;
                [this_pat_loc,~] = ismember(nan_indexes , temp_pattern,'rows');
                
                X = current_observations(~temp_pattern,this_pat_loc)';
                current_mu_temp = current_mu(~temp_pattern);
                current_sig_temp = current_sig(~temp_pattern,~temp_pattern);
                if ~isempty(~temp_pattern)
                    PDFs = mvnpdf(X, current_mu_temp(:)', current_sig_temp);
                else
                    PDFs = ones(size(X,1),1);
                end
                P_observ_cond_to_state_comp(sum(state_numbers(1:zee-1))+s,k,this_pat_loc) = P_all(zee,s ,k)*PDFs;
            end
            
        end
    end
end

P_observ_cond_to_state = squeeze( sum( P_observ_cond_to_state_comp , 2) );

if(size(P_observ_cond_to_state,2)==1&&size(P_observ_cond_to_state_comp,1)==1)
    P_observ_cond_to_state=P_observ_cond_to_state';
end



end
