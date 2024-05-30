
function [pi_0_lsim , coupling_tetha_IM , transition_matrices_IM ,  gmm_para_lsim ,  AIC, log_likelyhood, BIC, pi_steady] = em_lsim_svi( channels_observations , channel_num_states , num_gmm_component , max_itration , extra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVI inference
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



if(extra.plot==1)
    figure
    drawnow
    pause(0.01)
end

%Initializing all parameters
C = size(channels_observations , 1);

dim_observation = zeros(C , 1);

for zee = 1:C
    dim_observation(zee,1)  = size( channels_observations{zee,1} ,1);
end

coupling_tetha_IM = ones(C , C )/C;
transition_matrices_IM = ones( sum(channel_num_states) , sum(channel_num_states)  );
pi_0_lsim = zeros( sum(channel_num_states) , 1  );

% MU_all is N*C x Dimension
P_all = zeros( C , max(channel_num_states) , max(num_gmm_component)  );
mu_all = zeros(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component)   );
sigma_all =  ones(  sum(dim_observation)  , max(channel_num_states) , max(num_gmm_component) );

state_numbers_index =  ( [0;cumsum(channel_num_states(:))] );
temp_coupling = zeros(state_numbers_index(end)+1,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_trails = size(channels_observations,2);
length_observation(1,1)= length(channels_observations{1,1}(1,:));

for i=2:num_trails
    length_observation(i,1)= length(channels_observations{1,i}(1,:));
end

T_all = sum(length_observation);

all_observation =  zeros( sum(dim_observation) , T_all  );
dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );
trail_time_index =  ( [0;cumsum(length_observation(:))] );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for zee = 1:C
    
    Y = channels_observations{zee,1};
    
    for i=2:num_trails
        Y = cat(2, Y, channels_observations{zee,i});
    end
    
    
    all_observation(dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1),:)=Y;
    
    % Initialization of emmision probabilties
    disp('Initializing: ')
    init_gmm_numbers = [1,channel_num_states(zee)*num_gmm_component(zee)];
    max_Itr = 5;
    Replicates_Number = 20;
    
    temp_targets = ones(1,size(Y,2));
    
    warning('off')
    gmm_para  = fine_gmm_fitting( Y  , temp_targets , init_gmm_numbers , max_Itr , Replicates_Number );
    warning('on')
    
    temp_mu = gmm_para.x.mu;
    temp_sigma = gmm_para.x.Sigma;
    P = gmm_para.x.ComponentProportion;
    
    gmm_para=[];
    
    for i=1:channel_num_states(zee)
        
        
        mu_state = temp_mu(1 ,:)';
        mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , 1 ) = temp_mu(1 ,:)';
        sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , 1 ) = squeeze(temp_sigma(1,:,1));
        P_all(zee , i , 1) = P(1);
        
        P(1)=[];
        temp_mu(1 ,:)=[];
        temp_sigma(:,:,1)=[];
        
        
        for k=2:num_gmm_component(zee)
            temp_mu_gmm = mu_state(:,1);
            temp_distance = sum( abs(repmat(temp_mu_gmm,1,size(temp_mu,1)) - temp_mu').^2 );
            
            [~,ind_min]=min(temp_distance);
            
            
            mu_all(  dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i ,  k) = temp_mu(ind_min ,:)';
            sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k) = squeeze(temp_sigma(1,:,ind_min));
            P_all(zee , i , k)=P(ind_min);
            
            P(ind_min)=[];
            temp_mu(ind_min ,:)=[];
            temp_sigma(:,:,ind_min)=[];
            
        end
        
        P_all(zee , i , :)= P_all(zee , i , :)./( sum( P_all(zee , i , :) ,3)  );
        
    end
    
    
    % initial probabilities
    pi_0_lsim(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1) = ones(channel_num_states(zee) ,1  )/channel_num_states(zee);
    
    
    % Transition probabilities  M(c) x M(zee)
    transition_matrices_IM( : , state_numbers_index(zee)+1:state_numbers_index(zee+1) ) =  1/(-state_numbers_index(zee)+state_numbers_index(zee+1));
    
    
end




transitions_matrices_new = transition_matrices_IM;
coupling_tetha_new = coupling_tetha_IM ;
P_all_new = P_all;
mu_all_new = mu_all;
sigma_all_new =  sigma_all;



%  from here gpu begin
alpha_T_trails = zeros( sum(channel_num_states) , T_all  );
% Gamma_P0 = zeros( sum(State_Numbers) , Trails );

alpha_trails = zeros( sum(channel_num_states) , T_all );
beta_trails = zeros( sum(channel_num_states) , T_all );

b_c_ot_nc_trails = zeros( sum(channel_num_states) , T_all );
P_observ_cond_to_state_trails = zeros( sum(channel_num_states) , T_all  );
P_observ_cond_to_state_comp_trails = zeros( sum(channel_num_states) , max(num_gmm_component) , T_all  );

% Zeta_trails = zeros( sum(State_Numbers) , sum(State_Numbers) , T_all-Trails  );

gamma_trail_time_index = 1+(trail_time_index(1:end-1));
zeta_trail_time_index = 1:trail_time_index(end);
zeta_trail_time_index(trail_time_index(2:end))=[];
L_T_1 =  ( length(zeta_trail_time_index) );



P_O_model =  zeros( 1 , num_trails );


coupling_transition_all = zeros(state_numbers_index(end) , state_numbers_index(end)   );


mat_mult_one = zeros( state_numbers_index(end) , C  );

vec_mat_mult=[];
for zee =  1:C
    vec_mat_mult = cat(1,vec_mat_mult, [(zee-1)*state_numbers_index(end) + [state_numbers_index(zee)+1:state_numbers_index(zee+1)]]' );
end

mat_mult_one(vec_mat_mult)=1;


for itr = 1:max_itration
    
    
    for zee = 1:C
        
        zee_index = state_numbers_index(zee)+1:state_numbers_index(zee+1);
        coupling_transition_all(: , zee_index ) =  repmat(mat_mult_one*coupling_tetha_IM(:,zee),1,channel_num_states(zee)).*transition_matrices_IM(:,zee_index)  ;
        
    end
    
    
    
    
    
    sum3_zeta = 0;
    s_last = state_numbers_index(end);
    for tr = 1:num_trails
        
        tr_index = trail_time_index(tr)+1:trail_time_index(tr+1);
        tr_indexm1 = trail_time_index(tr)+1:trail_time_index(tr+1)-1;
        tr_indexp1 = trail_time_index(tr)+2:trail_time_index(tr+1);
        
        [P_O_model(tr) ,alpha_trails(:,tr_index) , beta_trails(:,tr_index)  , alpha_T_trails(:,tr_index)...
            , b_c_ot_nc_trails(:,tr_index) , P_observ_cond_to_state_trails(:,tr_indexv), P_observ_cond_to_state_comp_trails(:,:,tr_index) ] = ...
            forward_backward_chmm_kernel( pi_0_lsim , transition_matrices_IM , coupling_tetha_IM,  P_all , mu_all  , sigma_all ,   all_observation(:,tr_index) , state_numbers_index, channel_num_states ,  dimension_numbers_index , num_gmm_component, sigma_diag  );        

        zeta_trails = repmat( reshape(alpha_trails(:,tr_indexm1).*b_c_ot_nc_trails(:, tr_indexm1), s_last,1, length_observation(tr,1)-1) ,1, s_last,1)...
            .*repmat( reshape(beta_trails(:, tr_indexp1 ), 1 , s_last , length_observation(tr,1)-1), s_last,1,1);
        sum3_zeta = sum3_zeta + sum(zeta_trails,3);
        
    end
    
    sum3_zeta = sum3_zeta.* coupling_transition_all;
    
    %     zeta_trails = ( repmat( reshape(alpha_trails( : , zeta_trail_time_index ).*b_c_ot_nc_trails( : , zeta_trail_time_index ),state_numbers_index(end),1,L_T_1) ,1,state_numbers_index(end),1).*repmat( reshape(beta_trails(:,1+zeta_trail_time_index), 1 , state_numbers_index(end) , L_T_1),state_numbers_index(end),1,1) ) ;
    %     sum3_zeta = sum(zeta_trails,3).* coupling_transition_all;
    
    gamma_P0 = alpha_T_trails(: , gamma_trail_time_index);
    
    
    
    
    
    % update initial & transiotion probabilities
    
    pi_0_new = mean(gamma_P0,2);
    
    for zee = 1:C
        
        zee_index = state_numbers_index(zee)+1:state_numbers_index(zee+1);
        dim_zee = dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1);
        
        sum3_Zeta_zee = sum3_zeta(: , zee_index );
        sum32_Zeta_zee = sum(sum3_Zeta_zee,2);
        temp_Trans_zee  = sum3_Zeta_zee./repmat(sum32_Zeta_zee,1,channel_num_states(zee));
        temp_Trans_zee(isnan(temp_Trans_zee))=1/channel_num_states(zee);
        transitions_matrices_new(: , zee_index) = temp_Trans_zee;
        
        temp_coupling(2:end) = cumsum( sum32_Zeta_zee ,1 );
        coupling_tetha_new(:,zee)= (temp_coupling(state_numbers_index(2:end)+1)-temp_coupling(state_numbers_index(1:end-1)+1)) /L_T_1;
        
        % update emmision probabilities
        
        for i = 1:channel_num_states(zee)
            
            temp_gamma_obs = zeros(num_gmm_component(zee) , T_all );
            for k = 1:num_gmm_component(zee)
                
                temp_gamma_obs(k,:) = alpha_T_trails( zee_index(i),: ) .* (squeeze( P_observ_cond_to_state_comp_trails(zee_index(i),k,:) )'./P_observ_cond_to_state_trails( zee_index(i),:) );
                temp_gamma_obs(k,isnan(temp_gamma_obs(k,:)))=0;
                
                Y_update = all_observation(dim_zee,:);
                weight_update = repmat(temp_gamma_obs(k,:),dim_observation(zee,1),1);
                weight_update(isnan(Y_update))=0;
                Y_update(isnan(Y_update))=0;
                
                if sum(weight_update>0)==0
                    weight_update = weight_update+1;
                end
                mu_all_new( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k ) =  sum( weight_update.*Y_update , 2)./ sum( weight_update , 2);
                Recentering = Y_update - repmat(mu_all_new( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k ) , 1 , size(Y_update,2));
                
                sigma_all_new( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k ) =  0.001*var(all_observation(dim_zee,:),[],2,'omitnan') + ...
                    diag(  (weight_update.*Recentering)*Recentering')./sum( weight_update , 2) ;
                
                if(sum(temp_gamma_obs(:)<0)>0||sum(isnan(temp_gamma_obs(:)))>0)||sum(isnan(mu_all_new(:)))>0
                    disp('Error: reducing the state or gmm numbers')
                end
                
            end
            
            P_all_new(zee , i , 1:num_gmm_component(zee) ) = sum(temp_gamma_obs,2)/sum(temp_gamma_obs(:));
            
        end
        
        
    end
    
    pi_steady = sum( alpha_T_trails ,2)/ T_all;
    
    pi_0_lsim = pi_0_new;
    transition_matrices_IM = transitions_matrices_new;
    coupling_tetha_IM = coupling_tetha_new;
    P_all = P_all_new;
    mu_all = mu_all_new;
    sigma_all =  sigma_all_new;
    
    AIC(itr)= -2*sum(P_O_model) + 2* ( C*(C-1) +  sum(channel_num_states(:))*sum(channel_num_states(:)) - C*sum(channel_num_states(:)) + sum(channel_num_states(:).*num_gmm_component(:)*2.*dim_observation(:)) +  sum(channel_num_states(:).*(num_gmm_component(:)-1))  );
    BIC(itr)= -2*sum(P_O_model) + log(length(temp_targets))* ( C*(C-1) +  sum(channel_num_states(:))*sum(channel_num_states(:)) - C*sum(channel_num_states(:)) + sum(channel_num_states(:).*num_gmm_component(:)*2.*dim_observation(:)) +  sum(channel_num_states(:).*(num_gmm_component(:)-1))  );
    
    log_likelyhood(itr) = sum(P_O_model);
    
    %         disp(['Itr: ',num2str(itr),', Log likelyhood: ',num2str(sum(P_Y_Model)),', AIC: ',num2str(AIC(itr))])
    
    if(extra.check_convergence>0)
        if(itr>9)
            if(sum((diff(log_likelyhood(itr-5:itr)))<0.001)>=2)
                disp('EM algorithm converged')
                break;
            end
        end
    end
    
    if(itr>3)
        
        if(max_AIC<log_likelyhood(end))
            max_AIC = log_likelyhood(end);
            P_all_max = P_all;
            mu_all_max = mu_all;
            sigma_all_max =  sigma_all;
            transitions_matrices_max = transition_matrices_IM;
            coupling_tetha_max = coupling_tetha_IM;
            pi_0_max = pi_0_lsim;
            pi_steady_max = pi_steady;
        end
        
    else
        
        max_AIC = log_likelyhood(end);
        P_all_max = P_all;
        mu_all_max = mu_all;
        sigma_all_max =  sigma_all;
        transitions_matrices_max = transition_matrices_IM;
        coupling_tetha_max = coupling_tetha_IM;
        pi_0_max = pi_0_lsim;
        pi_steady_max = pi_steady;
        
    end
    
    if(rem(itr,10)==0&&extra.plot==1)
        hold on
        plot(log_likelyhood,'b')
        pause(0.01);
    end
    
end

if(extra.plot==1)
    hold on
    plot(log_likelyhood,'b')
    pause(0.01);
end


disp(['Itr: ',num2str(itr),', Log likelyhood: ',num2str(sum(P_O_model)),', AIC: ',num2str(AIC(end)),', BIC: ',num2str(BIC(end))])

%prepare output parameters

P_all = P_all_max;
mu_all = mu_all_max;
sigma_all =  sigma_all_max;

transitions_matrices_out = transitions_matrices_max;
pi_0_out = pi_0_max;
pi_steady_out = pi_steady_max;
coupling_tetha_IM = coupling_tetha_max;

clear Transitions_matrices PI_0 PI_Steady

gmm_para_lsim = cell(C,1);
transition_matrices_IM = cell(C,C);
pi_0_lsim = cell(C,1);
pi_steady = cell(C,1);

for zee = 1:C
    
    
    
    gmm_para=[];
    
    for i=1:channel_num_states(zee)
        
        for k=1:num_gmm_component(zee)
            
            gmm_para(i).P(k,1)=P_all(zee,i,k);
            gmm_para(i).sigma(k).x = sigma_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k);
            gmm_para(i).mu(k).x =  mu_all( dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1) ,  i , k);
            
        end
        
    end
    
    
    gmm_para_lsim{zee,1}.gmm_para = gmm_para;
    
    
    
    % initial probabilities
    pi_0_lsim{zee,1} =  pi_0_out(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1);
    pi_steady{zee,1} =  pi_steady_out(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1);
    
    for c=1:C
        % Transition probabilities  M(c) x M(zee)
        transition_matrices_IM{c,zee} =  transitions_matrices_out( state_numbers_index(c)+1:state_numbers_index(c+1) , state_numbers_index(zee)+1:state_numbers_index(zee+1) ) ;
    end
    
end


end



function gmm_para  = fine_gmm_fitting( inpute_features  , temp_targets , gmm_numbers , max_itr , replicates_number )




BaseBrices_Vectors_PDF = inpute_features;


options = statset('Display','off','MaxIter',max_itr);

class_Label = unique(temp_targets);

% minimum AIC value

for L = 1:length(class_Label)
    
    AIC_P=[];
    k=0;
    
    for num_Guassian = gmm_numbers
        k=k+1;
        
        flag1 = 0;
        counter = 0;
        while(flag1==0&&counter<10)
            counter = counter+1;
            try
                GMModel_P1 = fitgmdist( inpute_features(:,temp_targets==class_Label(L))' , num_Guassian ,'Start','randSample','Options', options,'Regularize',0.01*min(var(inpute_features(:,temp_targets==class_Label(L)),[],2,'omitnan')) , 'CovarianceType','diagonal','Replicates',replicates_number);
                GMModel_PP(k).x = GMModel_P1;
                AIC_P(k)=GMModel_P1.BIC;
                flag1=1;
                PostProbab = pdf( GMModel_P1 , BaseBrices_Vectors_PDF' )./( pdf( GMModel_P1,BaseBrices_Vectors_PDF' )+pdf( GMModel_P1 , BaseBrices_Vectors_PDF' ) );
            catch exception
                try
                    if(strfind(exception.message,'Ill-conditioned'))
                        disp('Requlized +1')
                        GMModel_P1 = fitgmdist( inpute_features(:,temp_targets==class_Label(L))' , num_Guassian ,'Start','randSample','Options', options,'Regularize',0.01*min(var(inpute_features(:,temp_targets==class_Label(L)),[],2,'omitnan')) , 'CovarianceType','diagonal','Replicates',replicates_number);
                        GMModel_PP(k).x = GMModel_P1;
                        AIC_P(k)=GMModel_P1.BIC;
                        PostProbab = pdf( GMModel_P1 , BaseBrices_Vectors_PDF' )./( pdf( GMModel_P1,BaseBrices_Vectors_PDF' )+pdf( GMModel_P1 , BaseBrices_Vectors_PDF' ) );
                        flag1=1;
                    end
                catch
                    
                    for i = 1:size(inpute_features,1)
                        
                        inpute_features(i,temp_targets==class_Label(L)) = inpute_features(i,temp_targets==class_Label(L))  + counter*0.05*mean(std(inpute_features(i,temp_targets==class_Label(L)),[],2))*randn(size(inpute_features(i,temp_targets==class_Label(L)) ,1) , size(inpute_features(i,temp_targets==class_Label(L)) ,2));
                        
                    end
                    disp('Not converge +1')
                    
                end
            end
        end
        
    end
    
    
    
    gmm_para(L).x =  GMModel_PP(2).x;
    gmm_para(L).BIC =  AIC_P;
    gmm_para(L).Num_Com = gmm_numbers(2);
    
    
end

end




function [P_O_model ,alpha , beta , alpha_T , b_c_ot_nc , P_observ_cond_to_state , P_observ_cond_to_state_comp]...
    = forward_backward_chmm_kernel( pi_0_chmm , transition_matrices_convex_comb , coupling_tetha_convex_comb, P_all , mu_all  , sigma_all ,  channel_time_series , state_numbers_index, channel_num_states ,  dimension_numbers_index , num_gmm_component  )

%  all_observation(:,Trail_Time_Index(tr)+1:Trail_Time_Index(tr+1))

% Y is observations which is a cell size Cx1  where Y{c,1} is the c'th subsystem observations

% C is the total number of subsystems

%Transition probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transitions_matrices is a cell array CxC wich cell (i,j) contain
% parameters P(v^j_t|v^i_t-1) (similar to Rabiner) (i->j)

%each cell (i,j) contain a matrix M(i)xM(j) which shows the transition
%probabilities subsystem i on subsystem j according to Rabiner notations  A=Transitions_matrices{i,j} is transition matrix A(n_j , n_i)= P(v^j_t = n_j | v^i_t = n_i ) & sum_n_j (A(i , j)) = 1
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
% CHMM_GMM_Param{j}.GMM_Param(n_j).P is  mixture wieght (probability) of the
% k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.GMM_Param(n_j).mu(k).x is  mean vector of the k'th GMM component related to hidden state (v^j_t = n_j)

% CHMM_GMM_Param{j}.GMM_Param(n_j).Sigma(k).x is  Covariance matrix of the k'th GMM component related to hidden state (v^j_t = n_j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C =  ( length( channel_num_states ));

%Computing Alpha^c_(t|t-1) for all subsystems


%Y is PxT, P is observation dimension & T is observations length

T =  ( size( channel_time_series , 2) );

P_ot_c_cond_past = zeros( C , T  );
b_c_ot_nc  = zeros( state_numbers_index(end) , T );

alpha = zeros( state_numbers_index(end) , T  );
alpha_T  = zeros( state_numbers_index(end) , T  );
beta  = zeros( state_numbers_index(end) , T );
a_alpha_b = zeros(state_numbers_index(end) , C , T-1 );




[P_observ_cond_to_state  ,   P_observ_cond_to_state_comp ]  = ...
    gmm_pdf_fast( P_all , mu_all  , sigma_all  , channel_time_series   ,  channel_num_states ,  dimension_numbers_index , num_gmm_component  );


% zee_index{zee,1} = cell{C,1};
mat_mult = zeros( state_numbers_index(end) , C  );

cuopling_repmat = zeros( state_numbers_index(end) , C  );
cuopling_repmat_beta = zeros( state_numbers_index(end) , C  );
weights_subsystems_opt = zeros( state_numbers_index(end) , C  );


coupling_transmat = zeros( state_numbers_index(end) , state_numbers_index(end)  );
diag_zero = ones( state_numbers_index(end) , state_numbers_index(end)  );


mat_mult_one = mat_mult;
vec_mat_mult = [];


for zee =  1:C
    
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    
    vec_mat_mult = cat(1,vec_mat_mult, ((zee-1)*state_numbers_index(end) + ( state_numbers_index(zee)+1:state_numbers_index(zee+1) ) )' );
    
    zee_coupling = coupling_tetha_convex_comb(:,zee)';
    cuopling_repmat( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    diag_zero( zee_index ,  zee_index )= 0 ;
    
    zee_coupling = coupling_tetha_convex_comb(zee,:);
    cuopling_repmat_beta( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    zee_coupling = coupling_tetha_convex_comb(zee,:)/sum(coupling_tetha_convex_comb(zee,:));
    weights_subsystems_opt( zee_index ,: ) = repmat( zee_coupling ,channel_num_states(zee),1);
    
    
    
    
end

weights_subsystems_opt(isnan(weights_subsystems_opt)) = 1/C;

for zee = 1:C
    zee_index = ( state_numbers_index(zee)+1:state_numbers_index(zee+1) );
    coupling_transmat( : ,  zee_index ) = repmat(cuopling_repmat_beta(:,zee) , 1 , length(zee_index) ).* transition_matrices_convex_comb(: , zee_index);
end

mat_mult_one(vec_mat_mult)=1;

alpha( :  , 1) = pi_0_chmm(:);

P_ot_c_cond_past(: , 1) = mat_mult_one'*(alpha(: , 1).* P_observ_cond_to_state(:,1));
b_c_ot_nc( : , 1) = P_observ_cond_to_state(:,1) ./ (mat_mult_one*P_ot_c_cond_past(: , 1));

%Compute Forward variables & Scaling coefficients
for t = (2:T)
    
    alpha_tm1_tm1 = alpha(: , t-1).*b_c_ot_nc(:,t-1);
    
    mat_mult(vec_mat_mult) = alpha_tm1_tm1;
    
    
    a_alpha_b( : , : ,t-1)= transition_matrices_convex_comb'*mat_mult ;%compute from t=1 to t=T-1
    a_alpha_b_coupled = a_alpha_b(  : , : ,t-1).*cuopling_repmat;
    alpha(: , t) = sum( a_alpha_b_coupled , 2);
    
    P_ot_c_cond_past(: , t) = mat_mult_one'*(alpha(: , t).* P_observ_cond_to_state(:,t));
    
    emmbed_P_ot_c_cond_past = (mat_mult_one*P_ot_c_cond_past(: , t));
    b_c_ot_nc( : , t) = P_observ_cond_to_state(:,t) ./ emmbed_P_ot_c_cond_past;
    
    
end



beta(: , T) = b_c_ot_nc(: , T);
alpha_T(: , T) =  alpha(: , T).*beta(: , T);

flag_nan = sum(sum(coupling_tetha_convex_comb,2)<=10^(-10))>0;


%Compute Backward variables
for t = (T-1:-1:1)
    
    
    mat_mult(vec_mat_mult) = beta(: , t+1);
    
    beta_C = transition_matrices_convex_comb *mat_mult ;%compute from t=1 to t=T-1
    
    %weighted average beta
    
    
    
    
    beta(: , t) = b_c_ot_nc( : , t).*(  sum( beta_C.*weights_subsystems_opt , 2) ) ;
    
    if(flag_nan)
        beta(isnan(beta(: , t)) , t)=  b_c_ot_nc( isnan(beta(: , t)) , t);
    end
    
    
    w_opt = sum(weights_subsystems_opt,2);
    w_opt2 = max(0,C/2-w_opt);
    
    %     beta(: , t) = (w_opt2.*b_c_ot_nc( : , t) + w_opt.*beta(: , t));
    
    alpha_T( : , t) = beta(: , t).*alpha(: , t);
    
    temp_sum = cumsum(alpha_T( : , t));
    normalizer = temp_sum(state_numbers_index(2:end))-[0;temp_sum(state_numbers_index(2:end-1))];
    normalizer_vec = mat_mult_one*normalizer(:);
    
    beta( : , t) = beta( : , t)./normalizer_vec;
    alpha_T( : , t) = alpha_T( : , t)./normalizer_vec;
    
    
end



P_O_model = sum( log(P_ot_c_cond_past(:)) );

if(isnan(P_O_model))
    error('failed to compute log likelyhood')
end


end

