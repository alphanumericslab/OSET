
function [pi_0_ehmm , coupling_tetha_ehmm ,  transition_ehmm  ,ehmm_gmm_para, index_matrix ] = chmm_cartesian_product( pi_0_chmm ,  transition_chmm  ,chmm_gmm_para  )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converting CHMM to equivalent HMM
% ehmm stand for equivalent HMM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



coupling_tetha_ehmm = 1;


C = size(pi_0_chmm , 1);

temp =1;
for c=1:C    
    temp = kron(temp , pi_0_chmm{c,1});
end

pi_0_ehmm{1,1} = temp;

if size(chmm_gmm_para{1,1}.gmm_para(1).sigma(1).x,2)==1
    sigma_diag = 1;
else
    sigma_diag = 0;
end

dim_observation = zeros(C , 1);
state_numbers = zeros(C , 1);
num_gmm_component =  zeros(C , 1);



for zee = 1:C
    
    dim_observation(zee,1)  = length( chmm_gmm_para{zee,1}.gmm_para(1).mu(1).x );
    state_numbers(zee,1)  = size( pi_0_chmm{zee,1} ,1);
    num_gmm_component(zee,1)  = length(chmm_gmm_para{zee,1}.gmm_para(1).P);    
    
end


A_cartesian = zeros( prod(state_numbers) , prod(state_numbers) );

for row_num = 1:size(A_cartesian,1)
    
    temp =1;    
    for c=1:C        
        temp = kron(temp , transition_chmm{c,1}(row_num,:));
    end
    
    A_cartesian(row_num , :) = temp;

end

transition_ehmm{1,1} = A_cartesian;


index_matrix = zeros(C , prod(state_numbers));
index_matrix_gmm = zeros(C , prod(num_gmm_component));

for c = 1:C
    
    tempIndex = state_numbers;
    tempIndex(1:c)=[];    
    tempRaw = kron(1:state_numbers(c)  , ones(1,prod(tempIndex)));    
    tempIndex = state_numbers;
    tempIndex(c:end)=[];
    index_matrix(c,:)=repmat(tempRaw , 1 , prod(tempIndex) );
        
    tempIndex = num_gmm_component;
    tempIndex(1:c)=[];    
    tempRaw = kron(1:num_gmm_component(c)  , ones(1,prod(tempIndex)));    
    tempIndex = num_gmm_component;
    tempIndex(c:end)=[];
    index_matrix_gmm(c,:)=repmat(tempRaw , 1 , prod(tempIndex) );

end



dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );


channel_num_states = prod(state_numbers);
num_gmm_component_hmm = prod(num_gmm_component);
P_all_hmm = ones( 1 , max(channel_num_states) , max(num_gmm_component_hmm)  );
mu_all_hmm = zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component_hmm)   );

if sigma_diag
    sigma_all_hmm =  zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component_hmm) );
else
    sigma_all_hmm =  zeros(  sum(dim_observation)  , sum(dim_observation) , max(channel_num_states) , max(num_gmm_component_hmm) );
end



for i = 1:channel_num_states
    this_set =index_matrix(:,i);
    for k=1:num_gmm_component_hmm
        this_gmm = index_matrix_gmm(:,k);
        for c=1:C
            mu_all_hmm(dimension_numbers_index(c)+1:dimension_numbers_index(c+1),i,k) = chmm_gmm_para{c, 1}.gmm_para(this_set(c)).mu(this_gmm(c)).x;
             if sigma_diag
                sigma_all_hmm( dimension_numbers_index(c)+1:dimension_numbers_index(c+1) ,  i , k) = chmm_gmm_para{c, 1}.gmm_para(this_set(c)).sigma(this_gmm(c)).x;
            else
               sigma_all_hmm( dimension_numbers_index(c)+1:dimension_numbers_index(c+1), dimension_numbers_index(c)+1:dimension_numbers_index(c+1), i , k) = chmm_gmm_para{c, 1}.gmm_para(this_set(c)).sigma(this_gmm(c)).x;
             end
             P_all_hmm(1,i,k) = P_all_hmm(1,i,k)* chmm_gmm_para{c, 1}.gmm_para(this_set(c)).P(this_gmm(c));
        end
    end
end




gmm_para=[];

for i=1:channel_num_states

    for k=1:num_gmm_component_hmm

        gmm_para(i).P(k,1)=P_all_hmm(1,i,k);
        if sigma_diag
            gmm_para(i).sigma(k).x = sigma_all_hmm( : ,  i , k)  ;
        else
            gmm_para(i).sigma(k).x = sigma_all_hmm( :,:, i , k)  ;
        end
        gmm_para(i).mu(k).x =  mu_all_hmm( : ,  i , k);
    end

end


ehmm_gmm_para{1,1}.gmm_para = gmm_para;





end







