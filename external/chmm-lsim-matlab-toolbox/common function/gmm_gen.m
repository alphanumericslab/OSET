function x = gmm_gen( gmm_param , state_number )

%GMM parameter include GMMParam().P , GMMParam().mu , GMMParam().Sigma 
% StateNumber is between 1 to number of hidden states

P = gmm_param(state_number).P;
n = Select_P( P );
x = diag(gmm_param(state_number).sigma(n).x)^0.5*randn(length(gmm_param(state_number).mu(n).x) , 1) + gmm_param(state_number).mu(n).x ;


