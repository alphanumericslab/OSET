

function [P_O_model ,alpha_out , beta_out , alpha_T_out , b_c_ot_nc_out , P_observ_cond_to_state_out , P_observ_cond_to_state_comp_out]...
    = forward_backward_chmm(pi_0_chmm, transition_matrices_chmm, chmm_gmm_para,  Y_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHMM forward-backward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Y is observations which is a cell size Cx1  where Y{c,1} is the c'th subsystem observations

% C is the total number of subsystems

%Transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transitions_matrices is a cell array CxC which cell (i,j) contain
% parameters P(v^j_t|v^i_t-1) (similar to Rabiner) (i->j)

%each cell (i,j) contain a matrix M(i)xM(j) which shows the transition
%probabilities subsystem i on subsystem j according to Rabiner notations  A=Transitions_matrices{i,j} is transition matrix A(n_j , n_i)= P(v^j_t = n_j | v^i_t = n_i ) & sum_n_j (A(i , j)) = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Coupling probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coupling matrice is a CxC  which  (c,zee) contain
% parameters Coupling_Tetha(c,zee)
%So sum(Coupling_Tetha)=ones(1,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initials hidden states probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PI_0{j} contain a vector with length M(j) which contain initials hidden
% states probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%emmision probabilities of CHMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHMM_GMM_Param{j}.gmm_para(n_j).P is  mixture wieght (probability) of the
% k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).mu(k).x is  mean vector of the k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.gmm_para(n_j).sigma(k).x is  Covariance matrix of the k'th GMM component related to hidden state (v^j_t = n_j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initializing all parameters
C = size(Y_in , 1);

dim_observation = zeros(C , 1);
channel_num_states = zeros(C , 1);
num_gmm_component =  zeros(C , 1);

for zee = 1:C
    dim_observation(zee,1)  = size( Y_in{zee,1} ,1);
    channel_num_states(zee,1)  = size( pi_0_chmm{zee,1} ,1);
    num_gmm_component(zee,1)  = length(chmm_gmm_para{zee,1}.gmm_para(1).P);
end

state_numbers_index =  ( [0;cumsum(channel_num_states(:))] );
dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );

% MU_all is N*C x Dimension
P_all = zeros( C , max(channel_num_states) , max(num_gmm_component)  );
mu_all = zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component)   );

zee_temp = find(dim_observation>1);

if isempty(zee_temp)
    zee_temp=1;
else
    zee_temp = zee_temp(1);
end

if size(chmm_gmm_para{zee_temp,1}.gmm_para(1).sigma(1).x,2)==1
    sigma_all =  ones(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component) );
    sigma_diag = 1;
else
    sigma_all =  ones(  sum(dim_observation)  , max(dim_observation) , max(channel_num_states) , max(num_gmm_component) );
    sigma_diag = 0;
end

for zee = 1:C

    gmm_para = chmm_gmm_para{zee,1}.gmm_para ;

    for i=1:channel_num_states(zee)
        for k=1:num_gmm_component(zee)
            P_all(zee,i,k) =  gmm_para(i).P(k,1);
            if sigma_diag
                sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) = gmm_para(i).sigma(k).x ;
            else
                sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1), 1:dim_observation(zee), i , k) = gmm_para(i).sigma(k).x ;
            end
            mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) =  gmm_para(i).mu(k).x;
        end
    end

end


%Computing Alpha^c_(t|t-1) for all subsystems

channel_time_series = zeros(sum(dim_observation) , size(Y_in{1,1},2));
%Y is PxT, P is observation dimension & T is observations length
for zee=1:C
    channel_time_series(dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1),:)= Y_in{zee,1};
end


[pi_0_ehmm, coupling_tetha_ehmm,  transition_ehmm, ehmm_gmm_para, index_matrix ] = chmm_cartesian_product( pi_0_chmm ,  transition_matrices_chmm  ,chmm_gmm_para  );
Y_ehmm{1,1} = channel_time_series;
%& exact solution
[P_O_model, alpha_eq , ~, alpha_T_eq] = forward_backward_lsim_svi( pi_0_ehmm , coupling_tetha_ehmm , transition_ehmm ,  ehmm_gmm_para  ,  Y_ehmm  );

T =  ( size( channel_time_series , 2) );
alpha = zeros( state_numbers_index(end) , T  );
alpha_T  = zeros( state_numbers_index(end) , T  );
b_c_ot_nc  = zeros( state_numbers_index(end) , T );% dummy output

tempT_exact = alpha_T_eq{1};
temp_exact = alpha_eq{1};

for c=1:C
    %extract exact solutions
    index_temp = index_matrix(c,:);
    for s=1:channel_num_states(c)
        alpha_T(state_numbers_index(c)+s, :) = sum(tempT_exact( index_temp==s , :),1);
        alpha(state_numbers_index(c)+s, :) = sum(temp_exact( index_temp==s , :),1);
    end    
end
beta = alpha_T./(10^-12+alpha);



if(nargout>1)

    if sigma_diag
        [P_observ_cond_to_state  ,   P_observ_cond_to_state_comp ]  = ...
            gmm_pdf_fast( P_all , mu_all  , sigma_all  , channel_time_series   ,  channel_num_states ,  dimension_numbers_index , num_gmm_component  );
    else
        [P_observ_cond_to_state  ,   P_observ_cond_to_state_comp ]  = ...
            gmm_pdf( P_all , mu_all  , sigma_all  , channel_time_series   ,  channel_num_states ,  dimension_numbers_index , num_gmm_component  );
    end

    alpha_out = cell(C,1);
    beta_out = cell(C,1);
    alpha_T_out = cell(C,1);
    b_c_ot_nc_out = cell(C,1);
    P_observ_cond_to_state_out = cell(C,1);
    P_observ_cond_to_state_comp_out = cell(C,1);

    for zee = 1:C

        alpha_out{zee,1} =  alpha(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        beta_out{zee,1} =  beta(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);

        alpha_out{zee,1} =  alpha(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        beta_out{zee,1} =  beta(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);

        alpha_T_out{zee,1} =  alpha_T(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        b_c_ot_nc_out{zee,1} =  b_c_ot_nc(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);

        P_observ_cond_to_state_out{zee,1} =  P_observ_cond_to_state(state_numbers_index(zee)+1:state_numbers_index(zee+1) , :);
        P_observ_cond_to_state_comp_out{zee,1} = P_observ_cond_to_state_comp(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1:num_gmm_component(zee) , :);

    end


end

end



