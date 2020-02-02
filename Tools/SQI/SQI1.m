function ind = SQI1(x, peaks)
% Channel selection based on Pseudo-periodicity
%
% Fahimeh Jamshidian Tehrani
% February 2020

% Comment out for individual using
% % %     Remove the mean
% % mn = mean(x,2)*ones(1,size(x,2));
% % x = x - mn;
% % 
% % %     Normalize the variance of x
% % x = x./var(x,0,2);

[CH, ~] = size(x);

for ch = 1: CH
    
    Peakloc{ch} = find(peaks{ch});
    i = 1;
    T = mean(diff(Peakloc{ch}));
    winL = round(T/2);
    while Peakloc{ch}(i)<winL
        RWAmatrix(i,1:2*winL) = [zeros(1,winL-Peakloc{ch}(i)) x(ch,1:Peakloc{ch}(i)+winL)];
        i = i +1;
    end
    for j = i : length(Peakloc{ch})-1
        if Peakloc{ch}(j)+winL <= length(x(ch, :))
            RWAmatrix(j, 1:2*winL) = x(ch,Peakloc{ch}(j)-winL+1 : Peakloc{ch}(j)+winL);
        else
            RWAmatrix(j, 1:2*winL) = [x(ch, Peakloc{ch}(j)-winL+1:end) zeros(1,2*winL-length(x(ch, Peakloc{ch}(j)-winL+1:end))) ];
        end
    end
    if 2*winL > length(x(ch, Peakloc{ch}(end)-winL+1:end))
        RWAmatrix(length(Peakloc{ch}), 1:2*winL) = [x(ch, Peakloc{ch}(end)-winL+1:end) zeros(1,2*winL-length(x(ch, Peakloc{ch}(end)-winL+1:end))) ];
    else
        RWAmatrix(length(Peakloc{ch}), 1:2*winL) = x(ch, Peakloc{ch}(end)-winL+1:Peakloc{ch}(end)+winL);
    end
    
    
    mk{ch} = RWAverage(RWAmatrix);
    
    ind(ch) = mean(mean((RWAmatrix - ones(size(RWAmatrix,1),1)* mk{ch}).^2));   
    
    
%     tt  = mean((RWAmatrix - ones(size(RWAmatrix,1),1)* mk{ch}).^2);
%     figure, 
%     fig = errorbar(mk{ch}, tt/2,'b'), hold on, plot(mk{ch},'r','linewidth',2);
%     grid on
%     xlim([0 length(tt)])
%     set(gca,'LooseInset',get(gca,'TightInset'))
% %     saveas(fig, ['.\newFigure\SQI1' num2str(ch) '.jpg']); saveas(fig, ['.\newFigure\SQI1' num2str(ch) '.eps']); saveas(fig, ['.\newFigure\SQI1' num2str(ch) '.fig']);
%     
%     figure,
%     fig = plot(mk{ch},'r','linewidth',2);
%     grid on
%     xlim([0 length(tt)])
%     set(gca,'LooseInset',get(gca,'TightInset'))
% %     saveas(fig, ['.\newFigure\pureSQI1' num2str(ch) '.jpg']); saveas(fig, ['.\newFigure\pureSQI1' num2str(ch) '.eps']); saveas(fig, ['.\newFigure\pureSQI1' num2str(ch) '.fig']);
    
    
clear RWAmatrix;
end