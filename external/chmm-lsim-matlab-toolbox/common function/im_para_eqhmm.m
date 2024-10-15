
function [pi_0_ehmm , coupling_tetha_ehmm ,  transition_ehmm  ,ehmm_gmm_para, index_matrix, pi_0_chmm ,  transition_chmm  , chmm_gmm_para] = im_para_eqhmm(pi_0_lsim, lsim_gmm_para, coupling_tetha_convex_comb, transition_matrices_convex_comb)

C = length(pi_0_lsim);
channel_num_states = zeros(1,C);

for c = 1:C
    channel_num_states(c)= length(pi_0_lsim{c});
end

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

%equal CHMM transition matrix of the convex combination model for channel zee
temp_var = zeros(C,1);
for zee = 1:C
    for i = 1:channel_num_states(zee)
        
        for j=1:size(index_matrix , 2)
            colomn_number_matrix_j = weight_state_column' *  ( index_matrix(:,j)-1) + 1;
            for c = 1:C
                temp_var(c,1) = transition_matrices_convex_comb{c,zee}(index_matrix(c,j),i);
            end
            transition_chmm{zee,1}( colomn_number_matrix_j , i  )= coupling_tetha_convex_comb(:,zee)' * temp_var(:);
        end
        
    end
end

pi_0_chmm = pi_0_lsim;
chmm_gmm_para = lsim_gmm_para;
% generating equivalent HMM parameters to perform exact inference
[ pi_0_ehmm , coupling_tetha_ehmm ,  transition_ehmm  ,ehmm_gmm_para  ] = chmm_cartesian_product( pi_0_chmm ,  transition_chmm  , chmm_gmm_para  );


