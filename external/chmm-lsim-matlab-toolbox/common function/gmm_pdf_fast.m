
function [P_observ_cond_to_state,P_observ_cond_to_state_comp] = gmm_pdf_fast( P_all , mu_all  , sigma_all  , observations_time_series   ,  state_numbers ,  dimension_numbers_index , num_gmm_component  )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function work for diogonal covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_Observ_cond_to_State( t , State ) is P(y_t | x_t= State) & P_Observ_cond_to_State is a matrix (K x T)
% z = mvnpdf(X,MU,SIGMA)

%GMM parameter include GMMParam(state).P , GMMParam(state).mu , GMMParam(state).Sigma

% StateNumber is between 1 to cardinality of hidden states

% global P_Observ_cond_to_State GMMParam StateNumbers Observations nan_index correct_index miss_index P_Observ_cond_to_State_Comp

indexvalid =[];
for zee =  1:length(state_numbers)    
    indexvalid = cat(1, indexvalid , ( (zee-1)*max(state_numbers)+1:(zee-1)*max(state_numbers)+state_numbers(zee) )');    
end

[~ , T_all]=   size(observations_time_series) ;

observations_time_series = repmat(observations_time_series , 1 , 1, max(state_numbers) , max(num_gmm_component));

mu_all = reshape(mu_all , size(mu_all,1) , 1 , size(mu_all , 2) , size(mu_all , 3) );
mu_all = repmat(mu_all , 1 , T_all , 1 ,1);


sigma_all = reshape(sigma_all , size(sigma_all,1) , 1 , size(sigma_all , 2) , size(sigma_all , 3) );
sigma_all = repmat(sigma_all , 1 , T_all , 1 ,1);

P_all = reshape(P_all , size(P_all,1) , 1 , size(P_all , 2) , size(P_all , 3) );
P_all = repmat(P_all , 1 , T_all , 1 ,1);


ind_nan = isnan(observations_time_series);
observations_time_series(ind_nan) = mu_all(ind_nan);
quadform = (observations_time_series - mu_all).^2 ./sigma_all;

logSqrtDetSigma = 0.5*log(sigma_all);

% R = sqrt(Sigma_all);
% quadform = ((Observations - MU_all)./R).^2 ;
% logSqrtDetSigma = log(R);

logSqrtDetSigma (ind_nan)=0;
% quadform = sum(xRinv.^2, 2);
%
% y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);

if(sum(diff(dimension_numbers_index)==1)==length(dimension_numbers_index)-1)
    
    d = 1;
    PDFs = P_all.*exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
    
else
    
    D = log(2*pi)/2;
    D = repmat(D , size(quadform,1) , size(quadform,2) , size(quadform,3) , size(quadform,4));
    D(ind_nan)  = 0;
    logvalue = cumsum(-0.5*quadform - logSqrtDetSigma - D);
    logvalue = cat(1,zeros(1, size(quadform,2) , size(quadform,3) , size(quadform,4)) , logvalue);
    PDFs = P_all.*exp( logvalue(dimension_numbers_index(2:end)+1 , :,:,:) - logvalue( dimension_numbers_index(1:end-1)+1 , :,:,:) );
    
%     logvalue = cat(1,zeros(1, size(quadform,2) , size(quadform,3) , size(quadform,4)) , logvalue);
%     temp_log = logvalue(dimension_numbers_index(2:end)+1 , :,:,:) - logvalue( dimension_numbers_index(1:end-1)+1 , :,:,:);
%     if sum(sum(sum(temp_log<-700,4)==state_numbers(:).*num_gmm_component(:),2)>0
%     for c = 1:length(state_numbers)
%     index_t = sum(temp_log(c,:,:)<-700,3)==state_numbers(c);
%     temp_log(c,:,1:state_numbers(c))
%     end
%     end

end

PDF_permuted = permute(PDFs , [3,1,4,2]);
P_observ_cond_to_state_comp = reshape(PDF_permuted , size(PDF_permuted,1)*size(PDF_permuted,2) , size(PDF_permuted,3) , size(PDF_permuted,4));
P_observ_cond_to_state_comp = P_observ_cond_to_state_comp(indexvalid,:,:);
P_observ_cond_to_state = squeeze( sum( P_observ_cond_to_state_comp , 2) );

if(size(P_observ_cond_to_state,2)==1&&size(P_observ_cond_to_state_comp,1)==1)
    P_observ_cond_to_state=P_observ_cond_to_state';
end



end
