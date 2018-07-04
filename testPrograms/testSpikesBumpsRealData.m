% test the Spikes and Bumps phenomenon in ICA using real fetal ECG data
% Reza Sameni, Copyright 2006
% Modified June 2018

clear;
close all;

load LTE07102005_40s.txt -mat;

data = LOUETTE_40s';clear LOUETTE_40s
fs = 1000;

% Keep the follwoing lines if you want to repeat the same results again
randn('seed', 1);
rand('seed', 1);
%//////////////////////////////////////////////////////////////////////////
% remove undesired data
I = [7,16,18,19,30,31,41,51,52,57,58,59];
data = RemoveChannels(data,I);
data = data(:,1:20*fs);
%//////////////////////////////////////////////////////////////////////////
data = data - LPFilter(data,1.5/fs);
% data = LPFilter(data,150/fs);
time = (0:length(data)-1)/fs;

%//////////////////////////////////////////////////////////////////////////
% [s, A, W] = fastica(data,'approach', 'symm','displayMode','off','lastEig',55);
[s, A, W] = fastica(data,'approach', 'symm','displayMode','off');
%//////////////////////////////////////////////////////////////////////////
PlotECG(s,10,'b', fs);

% s1 = s([51 1 46 12 39 9 18 16],:);
s1 = s([5 6 54 37 44 3 39 48],:);
%//////////////////////////////////////////////////////////////////////////
I = 1:20*fs;
L = 8;
k = 0;
% names = {'i','ii','iii','avr','avl','avf','v_1','v_2','v_3','v_4','v_5','v_6','v_x','v_y','v_z'};
names = {'mECG_1','mECG_2','fECG','n_1','n_2','b_1','b_2','b_3'};
for i = 1:size(s1,1)
    if(mod(i,L)==1)
        h = figure;
        k = k + 1;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(time(I),s1(i,I),'k');

    a = axis;
    a(1) = time(1);
    a(2) = 20;
    axis(a);
    ylabel(names{i},'FontSize',16);
    grid;
    set(gca,'Box','On','FontSize',16);
    if (mod(i,L)~=0)
        set(gca,'XTickLabel',[]);
    else
        xlabel('Time (s)','FontSize',16);
    end

    if(mod(i,L)==0)
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[.01 .01 2.5 8.5]);
        % % %         print('-dpng','-r600',['C:\Reza\ECGChannels_',num2str(k),'.png']);
        % % %         print('-deps','-r600',['C:\Reza\ECGChannels_',num2str(k),'.eps']);
    end
% % %     set(gcf,'Position',[146 157 776 928]);
end

%//////////////////////////////////////////////////////////////////////////
