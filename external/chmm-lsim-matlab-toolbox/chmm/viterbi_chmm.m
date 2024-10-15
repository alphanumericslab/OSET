


function [P_star_model , X_star]...
    = viterbi_chmm(  pi_0_chmm , coupling_tetha_convex_comb , transition_matrices_convex_comb ,  chmm_gmm_para  ,  Y_in  )

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
state_numbers = zeros(C , 1);
num_gmm_component =  zeros(C , 1);

for zee = 1:C

    dim_observation(zee,1)  = size( Y_in{zee,1} ,1);
    state_numbers(zee,1)  = size( pi_0_chmm{zee,1} ,1);
    num_gmm_component(zee,1)  = length(chmm_gmm_para{zee,1}.gmm_para(1).P);

end

state_numbers_index =  ( [0;cumsum(state_numbers(:))] );
dimension_numbers_index =  ( [0;cumsum(dim_observation(:))] );

% MU_all is N*C x Dimension
P_all = zeros( C , max(state_numbers) , max(num_gmm_component)  );
mu_all = zeros(  sum(dim_observation)  , max(state_numbers) , max(num_gmm_component)   );

if size(chmm_gmm_para{1,1}.gmm_para(1).sigma(1).x,2)==1
    sigma_all =  ones(  sum(dim_observation)  , max(state_numbers) , max(num_gmm_component) );
    sigma_diag = 1;
else
    sigma_all =  ones(  sum(dim_observation)  , max(dim_observation) , max(state_numbers) , max(num_gmm_component) );
    sigma_diag = 0;
end

transitions_matrices = ones( sum(state_numbers) , sum(state_numbers)  );
pi_0 = zeros( sum(state_numbers) , 1  );

for zee = 1:C

    gmm_para = chmm_gmm_para{zee,1}.gmm_para ;

    for i=1:state_numbers(zee)

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

    % initial probabilities
    pi_0(state_numbers_index(zee)+1:state_numbers_index(zee+1) , 1) = pi_0_chmm{zee,1};

    for c=1:C
        % Transition probabilities  M(c) x M(zee)
        transitions_matrices( state_numbers_index(c)+1:state_numbers_index(c+1) , state_numbers_index(zee)+1:state_numbers_index(zee+1) ) = transition_matrices_convex_comb{c,zee}  ;
    end

end




%Computing Alpha^c_(t|t-1) for all subsystems

Y_input = zeros(sum(dim_observation) , size(Y_in{1,1},2));
%Y is PxT, P is observation dimension & T is observations length
for zee=1:C
    Y_input(dimension_numbers_index(zee)+1:dimension_numbers_index(zee+1),:)= Y_in{zee,1};
end



%Computing Alpha^c_(t|t-1) for all subsystems


%Y is PxT, P is observation dimension & T is observations length

T =  ( size( Y_input , 2) );


%Viterbi variables & Best hidden path

X_star = zeros(C , T);

P_ot_c_cond_past = zeros( C , T  );

Alpha = zeros( state_numbers_index(end) , T  );
log_Phi  = zeros( state_numbers_index(end) , T  );
Si  = zeros( state_numbers_index(end) , T );

b_c_ot_nc  = zeros( state_numbers_index(end) , T );
a_alpha_b = zeros(state_numbers_index(end) , C , T-1 );


if sigma_diag
    P_Observ_cond_to_State  = gmm_pdf_fast( P_all , mu_all  , sigma_all, Y_input,  state_numbers,  dimension_numbers_index, num_gmm_component );
else
    P_Observ_cond_to_State  = gmm_pdf( P_all , mu_all  , sigma_all, Y_input,  state_numbers,  dimension_numbers_index, num_gmm_component );
end

% zee_index{zee,1} = cell{C,1};
mat_mult = zeros( state_numbers_index(end) , C  );
cuopling_repmat = zeros( state_numbers_index(end) , C  );
cuopling_repmat_zero = zeros( state_numbers_index(end) , C  );

mat_mult_one = mat_mult;
vec_Mat_mult = [];
digonal_Transmat_Transpose = zeros( state_numbers_index(end) , max(state_numbers)  );
Phi_repmat = -inf*zeros( max(state_numbers) ,  C );
indexPhi =[];

for zee =  1:C

    zee_index = state_numbers_index(zee)+1:state_numbers_index(zee+1) ;

    vec_Mat_mult = cat(1,vec_Mat_mult, ((zee-1)*state_numbers_index(end) + ( state_numbers_index(zee)+1:state_numbers_index(zee+1) ) )' );
    indexPhi = cat(1, indexPhi , ( (zee-1)*max(state_numbers)+1:(zee-1)*max(state_numbers)+state_numbers(zee) )' );
    Zee_Coupling = coupling_tetha_convex_comb(:,zee)';
    cuopling_repmat( zee_index ,: ) = repmat( Zee_Coupling ,state_numbers(zee),1);
    Zee_Coupling(zee)=0;
    cuopling_repmat_zero( zee_index ,: )=repmat(Zee_Coupling ,state_numbers(zee),1);
    digonal_Transmat_Transpose( zee_index ,  1:state_numbers(zee) ) = coupling_tetha_convex_comb( zee , zee)* transitions_matrices(zee_index , zee_index)';

end


mat_mult_one(vec_Mat_mult)=1;

%First step is computing forward variables
%Compute Forward variables & Scaling coefficients

Alpha( :  , 1) = pi_0(:);
P_ot_c_cond_past(: , 1) = mat_mult_one'*(Alpha(: , 1).* P_Observ_cond_to_State(:,1));
b_c_ot_nc( : , 1) = P_Observ_cond_to_State(:,1) ./ (mat_mult_one*P_ot_c_cond_past(: , 1));

%Compute Forward variables & Scaling coefficients
for t = (2:T)

    mat_mult(vec_Mat_mult)=(Alpha(: , t-1).*b_c_ot_nc(:,t-1));

    a_alpha_b( : , : ,t-1)= transitions_matrices'*mat_mult ;%compute from t=1 to t=T-1
    Alpha(: , t) = sum( a_alpha_b(  : , : ,t-1).*cuopling_repmat , 2);

    P_ot_c_cond_past(: , t) = mat_mult_one'*(Alpha(: , t).* P_Observ_cond_to_State(:,t));
    b_c_ot_nc( : , t) = P_Observ_cond_to_State(:,t) ./ (mat_mult_one*P_ot_c_cond_past(: , t));


end


%Initialization
log_Phi(: , 1) = log( b_c_ot_nc(: , 1).*pi_0(:) );
Si(: , 1) =  0;

%Compute Backward variables
for t = 2:T

    Phi_repmat(indexPhi) = log_Phi( : , t-1);

    Phi_repmat_out = repmat(Phi_repmat , 1 , 1 , max(state_numbers));
    Phi_repmat_out = permute(Phi_repmat_out , [1,3,2]);
    Phi_repmat_out = reshape(Phi_repmat_out , max(state_numbers) , 1 , C * max(state_numbers) );
    Phi_repmat_out = squeeze(Phi_repmat_out);

    Phi_repmat_out = Phi_repmat_out(: , indexPhi )';

    [max_res , arg_max] = max( log( digonal_Transmat_Transpose + repmat( sum(a_alpha_b(: ,  : ,t-1).*cuopling_repmat_zero , 2) , 1 ,max(state_numbers) ) ) +  Phi_repmat_out  ,[],2) ;%see relation 22

    log_Phi(: , t) = log( b_c_ot_nc(: , t) ) + max_res;
    Si(: , t) = arg_max;

end


P_star_model = zeros(C , 1);
Phi_repmat(indexPhi) = log_Phi( : , T);

% Phi_repmat_T = Phi_repmat';

[P_star_model(:,1) , X_star(:,T)] = max(  Phi_repmat  );


%Back-Tracking
for t = T-1:-1:1

    Phi_repmat(indexPhi) = Si( : , t+1);
    for zee=1:C
        X_star(zee,t) = Phi_repmat( X_star(zee,t+1) , zee )';
    end
end


P_star_model = sum(P_star_model);


end


