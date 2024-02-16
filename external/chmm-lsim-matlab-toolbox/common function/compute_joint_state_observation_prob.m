

function log_joint_prob = compute_joint_state_observation_prob( PI_0,  A, CHMM_GMM_Param, Y_in, X)



%Initializing all parameters
C = size(Y_in , 1);

Dim_Observation = zeros(C , 1);
State_Numbers = zeros(C , 1);
Num_GMM_Component =  zeros(C , 1);

for zee = 1:C
    
    Dim_Observation(zee,1)  = size( Y_in{zee,1} ,1);
    State_Numbers(zee,1)  = size( PI_0{zee,1} ,1);
    Num_GMM_Component(zee,1)  = length(CHMM_GMM_Param{zee,1}.GMM_Param(1).P);
        
end

State_Numbers_Index =  ( [0;cumsum(State_Numbers(:))] );
Dimension_Numbers_Index =  ( [0;cumsum(Dim_Observation(:))] );

% MU_all is N*C x Dimension
P_all = zeros( C , max(State_Numbers) , max(Num_GMM_Component)  );
MU_all = zeros(  sum(Dim_Observation)  , max(State_Numbers) , max(Num_GMM_Component)   );
Sigma_all =  ones(  sum(Dim_Observation)  , max(State_Numbers) , max(Num_GMM_Component) );


for zee = 1:C
    
    GMM_Param = CHMM_GMM_Param{zee,1}.GMM_Param ;
    
    for i=1:State_Numbers(zee)
        
        for k=1:Num_GMM_Component(zee)
            
            P_all(zee,i,k) =  GMM_Param(i).P(k,1);
            Sigma_all( Dimension_Numbers_Index(zee)+1:Dimension_Numbers_Index(zee+1) ,  i , k) = GMM_Param(i).Sigma(k).x ;
            MU_all( Dimension_Numbers_Index(zee)+1:Dimension_Numbers_Index(zee+1) ,  i , k) =  GMM_Param(i).mu(k).x;
            
        end
        
    end
    
    
    
    
end




%Computing Alpha^c_(t|t-1) for all subsystems

Y_input = zeros(sum(Dim_Observation) , size(Y_in{1,1},2));
%Y is PxT, P is observation dimension & T is observations length
for zee=1:C
    Y_input(Dimension_Numbers_Index(zee)+1:Dimension_Numbers_Index(zee+1),:)= Y_in{zee,1};
end


[P_Observ_cond_to_State   ]  = ...
    gmm_pdf_fast( P_all , MU_all  , Sigma_all  , Y_input   ,  State_Numbers ,  Dimension_Numbers_Index , Num_GMM_Component  );

P_Observ_cond_to_State_out = cell(C,1);
for zee = 1:C
    
    P_Observ_cond_to_State_out{zee,1} =  P_Observ_cond_to_State(State_Numbers_Index(zee)+1:State_Numbers_Index(zee+1) , :);
    
end



C = length(PI_0);
Channel_Num_States = zeros(1,C);

for c = 1:C
    Channel_Num_States(c) = size( PI_0{c,1} , 1);
end

WeightStateCal = cumprod( Channel_Num_States(C:-1:2) );
WeightStateCal=[1;WeightStateCal(:)];
WeightStateCal = WeightStateCal(C:-1:1);

log_joint_prob = 0;

for zee = 1:C
    %define all nescesary parameters including Forward & Backward variables & Scaling coefficients
    pi_0_c = PI_0{zee,1};
    log_joint_prob = log_joint_prob + log(pi_0_c(X(zee,1))) + log(P_Observ_cond_to_State_out{zee}(X(zee,1),1));
    
end


T = size(X,2);

for t=2:T
    
    ColomnNumber = WeightStateCal' *  ( X(: , t-1)-1) + 1;
    for zee=1:C
        
        log_joint_prob = log_joint_prob + log(  A{zee,1}( X(zee , t) , ColomnNumber) );
        log_joint_prob = log_joint_prob + log( P_Observ_cond_to_State_out{zee}(X(zee,t),t) );
        
    end
    
end



end



function [P_Observ_cond_to_State,P_Observ_cond_to_State_Comp] = gmm_pdf_fast( P_all , MU_all  , Sigma_all  , Observations   ,  State_Numbers ,  Dimension_Numbers_Index , Num_GMM_Component  )



% P_Observ_cond_to_State( t , State ) is P(y_t | x_t= State) & P_Observ_cond_to_State is a matrix (K x T)
% z = mvnpdf(X,MU,SIGMA)

%GMM parameter include GMMParam(state).P , GMMParam(state).mu , GMMParam(state).Sigma

% StateNumber is between 1 to cardinality of hidden states

% global P_Observ_cond_to_State GMMParam StateNumbers Observations nan_index correct_index miss_index P_Observ_cond_to_State_Comp

indexvalid =[];
for zee =  1:length(State_Numbers)    
    indexvalid = cat(1, indexvalid , ( (zee-1)*max(State_Numbers)+1:(zee-1)*max(State_Numbers)+State_Numbers(zee) )');    
end

[~ , T_all]=   size(Observations) ;

Observations = repmat(Observations , 1 , 1, max(State_Numbers) , max(Num_GMM_Component));

MU_all = reshape(MU_all , size(MU_all,1) , 1 , size(MU_all , 2) , size(MU_all , 3) );
MU_all = repmat(MU_all , 1 , T_all , 1 ,1);


Sigma_all = reshape(Sigma_all , size(Sigma_all,1) , 1 , size(Sigma_all , 2) , size(Sigma_all , 3) );
Sigma_all = repmat(Sigma_all , 1 , T_all , 1 ,1);

P_all = reshape(P_all , size(P_all,1) , 1 , size(P_all , 2) , size(P_all , 3) );
P_all = repmat(P_all , 1 , T_all , 1 ,1);


ind_nan = isnan(Observations);
Observations(ind_nan) = MU_all(ind_nan);
quadform = (Observations - MU_all).^2 ./Sigma_all;

logSqrtDetSigma = 0.5*log(Sigma_all);

% R = sqrt(Sigma_all);
% quadform = ((Observations - MU_all)./R).^2 ;
% logSqrtDetSigma = log(R);

logSqrtDetSigma (ind_nan)=0;
% quadform = sum(xRinv.^2, 2);
%
% y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);

if(sum(diff(Dimension_Numbers_Index)==1)==length(Dimension_Numbers_Index)-1)
    
    d = 1;
    PDFs = P_all.*exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
    
else
    
    D = log(2*pi)/2;
    D = repmat(D , size(quadform,1) , size(quadform,2) , size(quadform,3) , size(quadform,4));
    D(ind_nan)  = 0;
    logvalue = cumsum(-0.5*quadform - logSqrtDetSigma - D);
    logvalue = cat(1,zeros(1, size(quadform,2) , size(quadform,3) , size(quadform,4)) , logvalue);
    PDFs = P_all.*exp( logvalue(Dimension_Numbers_Index(2:end)+1 , :,:,:) - logvalue( Dimension_Numbers_Index(1:end-1)+1 , :,:,:) );
    
end

PDF_permuted = permute(PDFs , [3,1,4,2]);
P_Observ_cond_to_State_Comp = reshape(PDF_permuted , size(PDF_permuted,1)*size(PDF_permuted,2) , size(PDF_permuted,3) , size(PDF_permuted,4));
P_Observ_cond_to_State_Comp = P_Observ_cond_to_State_Comp(indexvalid,:,:);
P_Observ_cond_to_State = squeeze( sum( P_Observ_cond_to_State_Comp , 2) );

if(size(P_Observ_cond_to_State,2)==1&&size(P_Observ_cond_to_State_Comp,1)==1)
    P_Observ_cond_to_State=P_Observ_cond_to_State';
end



end








