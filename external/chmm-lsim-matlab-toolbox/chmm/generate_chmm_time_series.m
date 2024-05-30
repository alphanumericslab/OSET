

function [ channel_observation , channel_hidden_states ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para , transition_chmm , pi_0_chmm )


C = length(pi_0_chmm);

weight_state_column = cumprod( channel_num_states(C:-1:2) );
weight_state_column = [1;weight_state_column(:)];
weight_state_column = weight_state_column(C:-1:1);
channel_hidden_states=zeros(C , T);
channel_observation{C,1}=[];

for j =1:C
    
    channel_observation{j,1} = zeros(channel_dim_observ(j) , T);
    transition_chmm{j,1} =  transition_chmm{j,1}';
end


for j=1:C
    
    channel_hidden_states(j , 1) = select_p( pi_0_chmm{j,1} );
    channel_observation{j,1}(: , 1) = gmm_gen( chmm_gmm_para{j,1}.gmm_para, channel_hidden_states(j , 1) );
    
end


for t=2:T
    
    column_number = weight_state_column' *  ( channel_hidden_states(: , t-1)-1) + 1;
    for j=1:C
        channel_hidden_states(j , t) = select_p( transition_chmm{j,1}(: , column_number) );
        channel_observation{j,1}(: , t) = gmm_gen( chmm_gmm_para{j,1}.gmm_para , channel_hidden_states(j , t) );
    end
    
end


end


function x = select_p( pi_0 )

%This function select one state or guassian component in GMM
% pi_0 is vector such that sum is 1


temp = find(cumsum(pi_0) > rand(1,1));

x = temp(1);

end


function x = gmm_gen( gmm_param , state_number )

%GMM parameter include GMMParam().P , GMMParam().mu , GMMParam().Sigma
% StateNumber is between 1 to number of hidden states

P = gmm_param(state_number).P;
n = select_p( P );

if size(gmm_param(1).sigma(1).x,2)==1
    x = diag(gmm_param(state_number).sigma(n).x)^0.5*randn(length(gmm_param(state_number).mu(n).x) , 1) + gmm_param(state_number).mu(n).x ;
else
    x = (gmm_param(state_number).sigma(n).x)^0.5*randn(length(gmm_param(state_number).mu(n).x) , 1) + gmm_param(state_number).mu(n).x ;
end

end
