
function [scores_lsim, Label] = sim_calassification_lsim(C,T,train_number, test_number)

repeat_num = 1;

channel_num_states(1:C) = randi([2,6],1,C);
while(sum(channel_num_states)>25)
    channel_num_states(1:C) = randi([2,6],1,C);
end

channel_dim_observ(1:C) =  randi([1,5],1,C);
num_gmm_component(1:C) = randi([1,3],1,C);


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


clear transition_chmm transition_chmm_second

for zee = 1:C


    %pi_0 is the initial state probabilities for the channel zee

    temp_var = rand( channel_num_states(zee) , 1);
    temp_var = temp_var / sum(temp_var);
    pi_0_chmm{zee,1} = temp_var;

    temp_var = rand( channel_num_states(zee) , 1);
    temp_var = temp_var / sum(temp_var);
    pi_0_chmm_second{zee,1} = temp_var;



    temp = abs(randn( channel_num_states(zee) , prod(channel_num_states)));
    temp = temp ./ repmat( sum(temp) , channel_num_states(zee) , 1) ;
    transition_chmm{zee,1} = temp';

    temp = abs(randn( channel_num_states(zee) , prod(channel_num_states)));
    temp = temp ./ repmat( sum(temp) , channel_num_states(zee) , 1) ;
    transition_chmm_second{zee,1} = temp';

    % Gaussian mixture model initialization for each state in channel
    % zee
    for i=1:channel_num_states(zee)

        temp_var = rand(num_gmm_component(zee) , 1);
        temp_var = temp_var / sum(temp_var);
        chmm_gmm_para{zee,1}.gmm_para(i).P = temp_var;

        for k=1:num_gmm_component(zee)

            chmm_gmm_para{zee,1}.gmm_para(i).mu(k).x =  1*i+randn(channel_dim_observ(zee) , 1);
            chmm_gmm_para{zee,1}.gmm_para(i).sigma(k).x = 1+2*rand(channel_dim_observ(zee) , 1);

        end

        temp_var = rand(num_gmm_component(zee) , 1);
        temp_var = temp_var / sum(temp_var);
        chmm_gmm_para_second{zee,1}.gmm_para(i).P = temp_var;

        for k=1:num_gmm_component(zee)
            chmm_gmm_para_second{zee,1}.gmm_para(i).mu(k).x =  1*i+randn(channel_dim_observ(zee) , 1);
            chmm_gmm_para_second{zee,1}.gmm_para(i).sigma(k).x = 1+2*rand(channel_dim_observ(zee) , 1);
        end
    end
end


%generating observations time-series
clear train_obs_chmm_1  train_obs_hmm_1  train_obs_chmm_2  train_obs_hmm_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model 1
features_train_svm_1 = zeros(T*sum(channel_dim_observ),train_number);

for tr =1:train_number

    [ channels_observations , channel_hidden_states ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para , transition_chmm , pi_0_chmm );

    for i=1:C
        train_obs_chmm_1{i,tr} = channels_observations{i};
    end

    temp = cell2mat(channels_observations);
    train_obs_hmm_1{1,tr} = temp;
    features_train_svm_1(:,tr)=temp(:);

end


features_test_svm_1 = zeros(T*sum(channel_dim_observ),test_number);

for tr =1:test_number

    [ channels_observations , channel_hidden_states ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para , transition_chmm , pi_0_chmm );
    for i=1:C
        test_obs_chmm_1{i,tr} = channels_observations{i};
    end
    temp = cell2mat(channels_observations);
    test_obs_hmm_1{1,tr} = temp;
    features_test_svm_1(:,tr)=temp(:);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model 2
features_train_svm_2 = zeros(T*sum(channel_dim_observ),train_number);

for tr =1:train_number


    [ channels_observations_second , channel_hidden_states_second ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para_second , transition_chmm_second , pi_0_chmm_second );
    for i=1:C
        train_obs_chmm_2{i,tr} = channels_observations_second{i};
    end
    temp = cell2mat(channels_observations_second);
    train_obs_hmm_2{1,tr} = temp;
    features_train_svm_2(:,tr)=temp(:);

end



features_test_svm_2 = zeros(T*sum(channel_dim_observ),test_number);
for tr =1:test_number

    [ channels_observations_second , channel_hidden_states_second ] = generate_chmm_time_series( T , channel_num_states , channel_dim_observ , chmm_gmm_para_second , transition_chmm_second , pi_0_chmm_second );

    for i=1:C
        test_obs_chmm_2{i,tr} = channels_observations_second{i};
    end

    temp = cell2mat(channels_observations_second);
    test_obs_hmm_2{1,tr} = temp;
    features_test_svm_2(:,tr)=temp(:);

end

%% Training LSIM from observations for model 1

max_itration = 25;
extra.plot=0;
extra.check_convergence = 0;
extra.time_series = 0;


state_numbers_all = 2:6;
num_gmm_component_all = [ones(1,length(state_numbers_all)),2*ones(1,length(state_numbers_all)),3*ones(1,length(state_numbers_all))];
state_numbers_all = [state_numbers_all,state_numbers_all,state_numbers_all];

max_reinit = 3;
max_s = length(state_numbers_all);
clear AIC_par pi_0_chmm_par coupling_tetha_convex_comb_par transition_matrices_convex_par chmm_gmm_para_par

AIC_par = zeros(max_reinit,max_s);
pi_0_chmm_par{max_reinit,max_s}=[];
coupling_tetha_convex_comb_par{max_reinit,max_s}=[];
transition_matrices_convex_par{max_reinit,max_s}=[];
chmm_gmm_para_par{max_reinit,max_s}=[];


for s = 1:length(state_numbers_all)

    Channel_Num_States_t =state_numbers_all(s)*ones(1,C);
    num_gmm_component_temp(1:C) = num_gmm_component_all(s);

    for reinit_num = 1:max_reinit


        [pi_0_chmm_temp , coupling_tetha_convex_comb_temp , transition_matrices_convex_comb_temp ,  chmm_gmm_para_temp ,  AIC, log_likelyhood , BIC ,pi_steady] = ...
            em_lsim( train_obs_chmm_1 , Channel_Num_States_t , num_gmm_component_temp , max_itration , extra);

        AIC_par(reinit_num,s) = AIC(end);
        pi_0_chmm_par{reinit_num,s} = pi_0_chmm_temp;
        coupling_tetha_convex_comb_par{reinit_num,s} = coupling_tetha_convex_comb_temp;
        transition_matrices_convex_par{reinit_num,s} = transition_matrices_convex_comb_temp;
        chmm_gmm_para_par{reinit_num,s} = chmm_gmm_para_temp;

    end


end

[ ~ , index_min] = min(AIC_par(:));
[I_row, I_col] = ind2sub(size(AIC_par),index_min);

pi_0_1 = pi_0_chmm_par{I_row,I_col};
coupling_tetha_1 = coupling_tetha_convex_comb_par{I_row,I_col};
transitions_matrices_1 = transition_matrices_convex_par{I_row,I_col};
gmm_para_1 = chmm_gmm_para_par{I_row,I_col};

% Training CHMM from observations for model 2

for s = 1:length(state_numbers_all)

    Channel_Num_States_t =state_numbers_all(s)*ones(1,C);
    num_gmm_component_temp(1:C) = num_gmm_component_all(s);

    for reinit_num = 1:max_reinit

        [pi_0_chmm_temp , coupling_tetha_convex_comb_temp , transition_matrices_convex_comb_temp ,  chmm_gmm_para_temp ,  AIC, log_likelyhood , BIC ,pi_steady] = ...
            em_lsim( train_obs_chmm_2 , Channel_Num_States_t , num_gmm_component_temp , max_itration , extra);

        AIC_par(reinit_num,s) = AIC(end);
        pi_0_chmm_par{reinit_num,s} = pi_0_chmm_temp;
        coupling_tetha_convex_comb_par{reinit_num,s} = coupling_tetha_convex_comb_temp;
        transition_matrices_convex_par{reinit_num,s} = transition_matrices_convex_comb_temp;
        chmm_gmm_para_par{reinit_num,s} = chmm_gmm_para_temp;

    end


end

[ ~ , index_min] = min(AIC_par(:));
[I_row, I_col] = ind2sub(size(AIC_par),index_min);

pi_0_2 = pi_0_chmm_par{I_row,I_col};
Coupling_Tetha_2 = coupling_tetha_convex_comb_par{I_row,I_col};
transitions_matrices_2 = transition_matrices_convex_par{I_row,I_col};
gmm_para_2 = chmm_gmm_para_par{I_row,I_col};

% Forward-Backward Test

for tr =1:test_number

    for c = 1:C
        Y_test{c,1}=test_obs_chmm_1{c,tr};
    end

    P_Y_Model_1  = forward_backward_lsim( pi_0_1 , coupling_tetha_1  , transitions_matrices_1 ,  gmm_para_1 , Y_test );
    P_Y_Model_2 = forward_backward_lsim( pi_0_2 , Coupling_Tetha_2 , transitions_matrices_2   , gmm_para_2 , Y_test );

    F1(tr) = P_Y_Model_1-P_Y_Model_2;
end


for tr =1:test_number

    for c = 1:C
        Y_test{c,1}=test_obs_chmm_2{c,tr};
    end

    P_Y_Model_1  = forward_backward_lsim( pi_0_1 , coupling_tetha_1  , transitions_matrices_1 ,  gmm_para_1 , Y_test );
    P_Y_Model_2 = forward_backward_lsim(  pi_0_2 , Coupling_Tetha_2 , transitions_matrices_2   , gmm_para_2 , Y_test );

    F2(tr) = P_Y_Model_1-P_Y_Model_2;
end


scores_lsim(repeat_num,:) = [F1,F2];
Label = [ones(1,test_number),zeros(1,test_number)];


