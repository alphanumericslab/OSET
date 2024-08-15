function morphology_aligner(stacked_events, reference_event, delta)

% test VCG loop alignment
% UNDER TEST
% Reza Sameni, 2005

tau = -delta:delta;
num_channels = size(stacked_events, 1);
num_events = size(stacked_events, 2);
event_len = size(stacked_events, 3);

%//////////////////////////////////////////////////////////////////////////
er = zeros(1, length(tau));
for k = 1 : num_events
    Z = squeeze(stacked_events(:, k, :));
    for d = 1 : length(tau)
        J = [zeros(delta+tau(d)-1,N) ; eye(N) ; zeros(delta-tau(d)+1,N)];
        C = Z*J'*reference_event';
        [U,~,V] = svd(C);
        Q = U*V';
        alpha = trace(Z'*Q*reference_event*J)/trace(J'*(reference_event'*reference_event)*J);
        X = Z - alpha*Q*reference_event*J;
        er(tau+delta) = sqrt(sum(diag(X'*X)));
    end
    [mn,I] = min(er);
    tauhat(k) = I-delta;
    J = [zeros(delta+tauhat(k)-1,N) ; eye(N) ; zeros(delta-tauhat(k)+1,N)];
    C = Z*J'*reference_event';
    [U,S,V] = svd(C);
    Qhat(:,:,k) = U*V';
    Alphahat(k) = trace(Z'*Qhat(:,:,k)*reference_event*J)/trace(J'*reference_event'*reference_event*J);
end
%//////////////////////////////////////////////////////////////////////////
for k = 1:cntr
    Z(:,:,k) = Alphahat(k)*Qhat(:,:,k)*reference_event*[zeros(delta+tauhat(k)-1,N) ; eye(N) ; zeros(delta-tauhat(k)+1,N)];
end

figure
plot3(squeeze(Z(1,:,:)),squeeze(Z(2,:,:)),squeeze(Z(3,:,:)));grid;title('after alignment');
% figure
% plot3(squeeze(Z(1,:,:)),squeeze(Z(2,:,:)),squeeze(Z(3,:,:)));grid

%plot(squeeze(Z(3,:,3:end-3)));grid

figure
for ch = 1:3
    subplot(1,3,ch)
    plot(std(Z(ch,:,3:end-3),[],3))
    hold on;
    plot(std(AllECG(:,ch,3:end-3),[],3),'r')
    grid
end

