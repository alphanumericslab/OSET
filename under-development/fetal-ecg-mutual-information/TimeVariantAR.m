function [ARCoefs , ARCoefs2 , P]= TimeVariantAR(data,order,x0,q,R,p0,alpha,beta,gamma,Vmean);

N = length(data);
A = 1;
Q = q*eye(order);
Wmean = zeros(order,1);
ARCoefs = zeros(order,N);
Xminus = x0;    %zeros(order,1);
Pminus = p0*eye(order);
P = zeros(N,1);

% Pbar = zeros(order,order,N);
% Phat = zeros(order,order,N);
% Xhat = zeros(order,N);

% Filtering
for i = 1:N,
%     Pbar(:,:,i) = Pminus;
%     Xbar(:,i) = Xminus;

    % % %     %if(mod(i,1)==0),
    % % %     %if(abs(data(i))>=.1),

    %%%    if(norm(Pminus-Q,'fro')>q),

    if(i<order+1)
        H = [-data(i-1:-1:1)' zeros(1,order-i+1)];
    else
        H = -data(i-1:-1:i-order)';
    end

    Yminus = H * Xminus + Vmean;
    inov = data(i)-Yminus;
    
    %R = beta*R + (1-beta)*inov.^2;
    %Q = gamma*R*eye(order);

    K = Pminus * H'/(H * Pminus * H' + alpha*R);
    Pplus = ( (eye(order) - K * H) * Pminus * (eye(order) - K * H)' + K * R * K' )/alpha;
    Xplus = Xminus + K*(inov);                                % A posteriori state estimate

    % % %     else
    % % %         Pplus = Pminus;
    % % %         Xplus = Xminus;
    % % %     end

    Xminus = A * Xplus + Wmean;                              % State update
    Pminus = A * Pplus * A' + Q;

    ARCoefs(:,i) = Xplus;
    P(i) = max(diag(Pplus));

%     Phat(:,:,i) = Pplus;
%     Xhat(:,i) = Xplus;
end

% Smoothing
% PSmoothed = zeros(size(Phat));
% X = zeros(size(Xhat));
% PSmoothed(:,:,N) = Phat(:,:,N);
% X(:,N) = Xhat(:,N);
% for k = N-1:-1:1,
%     S = Phat(:,:,k) * A' /Pbar(:,:,k+1);
%     X(:,k) = Xhat(:,k) + S * (X(:,k+1) - Xbar(:,k+1)); % ?
%     PSmoothed(:,:,k) = Phat(:,:,k) - S * (Pbar(:,:,k+1) - PSmoothed(:,:,k+1)) * S';
% end
ARCoefs2 = ARCoefs;
