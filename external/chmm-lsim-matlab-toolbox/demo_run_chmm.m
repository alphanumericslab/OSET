clc
clear
close all

addpath(genpath(pwd))

%%
% This demo generates one chmm model observation, then compare our approximate inference with the exact inference

%C is the number of channels in CHMM
C = 3;
T = 2000; % T is the number of time samples


channel_num_states(1:C) = randi([2,6],1,C);
while(sum(channel_num_states)>25)
    channel_num_states(1:C) = randi([2,6],1,C);
end

channel_num_states_cum = cumsum(channel_num_states);

channel_dim_observ(1:C) =  randi([1,5],1,C);
num_gmm_component(1:C) = 3;

% Testing forward-backward algorithm

%indexing of Cartesian product
index_matrix = zeros(C , prod(channel_num_states));

for i = 1:C
    
    tempIndex = channel_num_states;
    tempIndex(1:i)=[];
    tempRaw = kron(1:channel_num_states(i)  , ones(1,prod(tempIndex)));
    tempIndex = channel_num_states;
    tempIndex(i:end)=[];
    index_matrix(i,:)=repmat(tempRaw , 1 , prod(tempIndex) );
    
end

weight_state_column = cumprod( channel_num_states(C:-1:2) );
weight_state_column=[1;weight_state_column(:)];
weight_state_column = weight_state_column(C:-1:1);



% CHMM random parameters initialization
coupling_tetha_convex_comb = abs(rand(C,C));
coupling_tetha_convex_comb = coupling_tetha_convex_comb./repmat( sum(coupling_tetha_convex_comb) , C , 1 ) ;

sparsity_para = randi(max(1,C))-1;
%     sparsity_para = 0;
coupling_tetha_convex_comb(1:sparsity_para,:) = 0.01*coupling_tetha_convex_comb(1:sparsity_para,:);
coupling_tetha_convex_comb = coupling_tetha_convex_comb./repmat( sum(coupling_tetha_convex_comb) , C , 1 ) ;


for zee = 1:C
    
    
    %pi_0 is the initial state probabilities for the channel zee
    temp_var = rand( channel_num_states(zee) , 1);
    temp_var = temp_var / sum(temp_var);
    pi_0_chmm{zee,1} = temp_var;
    
    %convex combination conditional transition matrices for the channel zee
    for c = 1:C
        
        temp_var = 0.001+abs(rand( channel_num_states(c) , channel_num_states(zee) ) );
        temp_var = temp_var.^(C);
        temp_var = temp_var ./ repmat( sum(temp_var,2) , 1 , channel_num_states(zee) ) ;
        temp_var =  exp(-3*coupling_tetha_convex_comb(c,zee))+temp_var;
        temp_var = temp_var ./ repmat( sum(temp_var,2) , 1 , channel_num_states(zee) ) ;
        transition_matrices_convex_comb{c,zee} = temp_var;
        
    end
    
    
    %equal CHMM transition matrix of the convex combination model for channel zee
    temp_var = zeros(C,1);
    
    for i = 1:channel_num_states(zee)
        
        for j=1:size(index_matrix , 2)
            colomn_number_matrix = weight_state_column' *  ( index_matrix(:,j)-1) + 1;
            for c = 1:C
                temp_var(c,1) = transition_matrices_convex_comb{c,zee}(index_matrix(c,j),i);
            end
            transition_matrices_chmm{zee,1}( colomn_number_matrix , i  )= coupling_tetha_convex_comb(:,zee)' * temp_var(:);
        end
        
    end
    
    
    % Gaussian mixture model initialization for each state in channel
    % zee
    for i=1:channel_num_states(zee)
        
        temp_var = rand(num_gmm_component(zee) , 1);
        temp_var = temp_var / sum(temp_var);
        chmm_gmm_para{zee,1}.gmm_para(i).P = temp_var;
        
        for k=1:num_gmm_component(zee)
            mu_temp = 1*i+randn(channel_dim_observ(zee) , 1);
            sig_temp = diag(1+2*rand(channel_dim_observ(zee),1)) + 0.5*ones(channel_dim_observ(zee),channel_dim_observ(zee));
            chmm_gmm_para{zee,1}.gmm_para(i).mu(k).x =  mu_temp;
            chmm_gmm_para{zee,1}.gmm_para(i).sigma(k).x = sig_temp;
        end
        
    end
    
end


% generating equivalent HMM parameters to perform exact inference
[ channels_observations , channel_hidden_states ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para , transition_matrices_chmm , pi_0_chmm );

%%


%2 combining svi-sed
[P_O_model ,alpha , beta , alpha_T ]...
    = forward_backward_lsim( pi_0_chmm , coupling_tetha_convex_comb , transition_matrices_convex_comb,  chmm_gmm_para  ,  channels_observations  );

% % %3 our svi
% % [P_O_model_svi ,alpha_svi , beta_svi , alpha_T_svi ]...
% %     = forward_backward_lsim_svi( pi_0_chmm , coupling_tetha_convex_comb , transition_matrices_convex_comb,  chmm_gmm_para  ,  channels_observations  );
% % 
% % %4 proposed sed
% % [P_O_model_sed ,alpha_sed , beta_sed , alpha_T_sed ]...
% %     = general_forward_backward_lsim( pi_0_chmm , coupling_tetha_convex_comb , transition_matrices_convex_comb,  chmm_gmm_para  ,  channels_observations  );

%%

extra.sigma_diag = 0;
extra.plot=1;
extra.check_convergence =0;
extra.time_series = 0;
max_itration = 100;

% channels_observations{2,1}(1,[1,5,20,100,600])=nan;

num_gmm_component(:)=2;
channel_num_states(:) = 5;

[pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim ,  AIC, log_likelyhood , BIC ,pi_steady] = ...
    em_lsim( channels_observations , channel_num_states , num_gmm_component , max_itration , extra);

[P_O_lsim ,alpha_out , beta_out , alpha_T_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
    = forward_backward_lsim( pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim  ,  channels_observations  );

[pi_0_chmm,  transition_matrices_chmm, chmm_gmm_para,  log_likelyhood, AIC1, BIC1] = em_chmm( channels_observations , channel_num_states , num_gmm_component , max_itration , extra);

[P_O_chmm ,alpha_out , beta_out , alpha_T_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
    = forward_backward_chmm(pi_0_chmm, transition_matrices_chmm, chmm_gmm_para,  channels_observations);

