clear;
close all;

fs = 1000;
data = load('s0010_re.txt');
% data = load('s0014lre.txt');


time = data(:,1)';
data = data(:,2:end)';

data = data - LPFilter(data,1.5/fs);
data = LPFilter(data,150/fs);

W = jadeR(data);
s = W*data;

[V,D] = eig(cov(data'));
s2 = V'*data;
s2 = s2(end:-1:1,:);

peaks = PeakDetection(data(13,:),1/fs);
[s3,W,A] = PiCA(data,peaks);

s4 = PiCA(s2(1:11,:),peaks);


%//////////////////////////////////////////////////////////////////////////
I = 1:5*fs;
L = 4;
k = 0;
names = {'i','ii','iii','avr','avl','avf','v_1','v_2','v_3','v_4','v_5','v_6','v_x','v_y','v_z'};
for i = 1:size(s,1)
    if(mod(i,L)==1)
        h = figure;
        k = k + 1;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(time(I),data(i,I),'k');
    
    a = axis;
    a(1) = time(1);
    a(2) = time(I(end));
    axis(a);
    ylabel(names{i},'FontSize',16);
    grid;
    set(gca,'Box','On','FontSize',16);
    if(mod(i,L)==0 || i==size(s,1))
        xlabel('Time (s)','FontSize',16);
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[.01 .01 3.5 8.5])
        %         print('-dpng','-r600',['C:\Reza\ECGChannels_',num2str(k),'.png']);
        %         print('-deps','-r600',['C:\Reza\ECGChannels_',num2str(k),'.eps']);
    else
        set(gca,'XTickLabel',[]);
    end
end

%//////////////////////////////////////////////////////////////////////////
% ICA

I = 1:5*fs;
L = 4;
k = 0;
for i = 1:size(s,1)
    if(mod(i,L)==1)
        h = figure;
        k = k + 1;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(time(I),s(i,I),'k');
    
    a = axis;
    a(1) = time(1);
    a(2) = time(I(end));
    axis(a);
    ylabel(['IC_{',num2str(i),'}'],'FontSize',16);
    grid;
    set(gca,'Box','On','FontSize',16);
    if(mod(i,L)==0 || i==size(s,1))
        xlabel('Time (s)','FontSize',16);
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[.01 .01 3.5 8.5])
        %         print('-dpng','-r600',['C:\Reza\ICChannels_',num2str(k),'.png']);
        %         print('-deps','-r600',['C:\Reza\ICChannels_',num2str(k),'.eps']);
    else
        set(gca,'XTickLabel',[]);
    end
end

%//////////////////////////////////////////////////////////////////////////
% PCA
I = 1:5*fs;
L = 4;
k = 0;
for i = 1:size(s,1)
    if(mod(i,L)==1)
        h = figure;
        k = k + 1;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(time(I),s2(i,I),'k');
    
    a = axis;
    a(1) = time(1);
    a(2) = time(I(end));
    axis(a);
    ylabel(['PC_{',num2str(i),'}'],'FontSize',16);
    grid;
    set(gca,'Box','On','FontSize',16);
    if(mod(i,L)==0 || i==size(s,1))
        xlabel('Time (s)','FontSize',16);
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[.01 .01 3.5 8.5])
        %         print('-dpng','-r600',['C:\Reza\PCChannels_',num2str(k),'.png']);
        %         print('-deps','-r600',['C:\Reza\PCChannels_',num2str(k),'.eps']);
    else
        set(gca,'XTickLabel',[]);
    end
end

%//////////////////////////////////////////////////////////////////////////
% PiCA
I = 1:5*fs;
L = 4;
k = 0;
for i = 1:size(s,1)
    if(mod(i,L)==1)
        h = figure;
        k = k + 1;
    end
    subplot(L,1,mod(i-1,L)+1);
    plot(time(I),s3(i,I),'k');
    
    a = axis;
    a(1) = time(1);
    a(2) = time(I(end));
    axis(a);
    ylabel(['\pic_{',num2str(i),'}'],'FontSize',16);
    grid;
    set(gca,'Box','On','FontSize',16);
    
    if(mod(i,L)==0 || i==size(s,1))
        xlabel('Time (s)','FontSize',16);
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[.01 .01 3.5 8.5])
        %         print('-dpng','-r600',['C:\Reza\PiCChannels_',num2str(k),'.png']);
        %         print('-deps','-r600',['C:\Reza\PiCChannels_',num2str(k),'.eps']);
    else
        set(gca,'XTickLabel',[]);
    end
end

% % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % % PiCA 2
% % % % % I = 1:5*fs;
% % % % % L = 4;
% % % % % k = 0;
% % % % % for i = 1:size(s4,1),
% % % % %     if(mod(i,L)==1)
% % % % %         h = figure;
% % % % %         k = k + 1;
% % % % %     end
% % % % %     subplot(L,1,mod(i-1,L)+1);
% % % % %     plot(time(I),s4(i,I),'k');
% % % % %
% % % % %     a = axis;
% % % % %     a(1) = time(1);
% % % % %     a(2) = time(I(end));
% % % % %     axis(a);
% % % % %     ylabel(['\pic_{',num2str(i),'}'],'FontSize',16);
% % % % %     grid;
% % % % %     set(gca,'Box','On','FontSize',16);
% % % % %
% % % % %     if(mod(i,L)==0 || i==size(s,1))
% % % % %         xlabel('Time (s)','FontSize',16);
% % % % %         set(h,'PaperUnits','inches');
% % % % %         set(h,'PaperPosition',[.01 .01 3.5 8.5])
% % % % % %         print('-dpng','-r600',['C:\Reza\PiCChannels_',num2str(k),'.png']);
% % % % % %         print('-deps','-r600',['C:\Reza\PiCChannels_',num2str(k),'.eps']);
% % % % %     else
% % % % %         set(gca,'XTickLabel',[]);
% % % % %     end
% % % % % end

%//////////////////////////////////////////////////////////////////////////

p = [13 14 15];
h = figure;
plot3(data(p(1),I),data(p(2),I),data(p(3),I),'b.');
set(gca,'Box','On','FontSize',14);
xlabel(names{p(1)},'FontSize',14);
ylabel(names{p(2)},'FontSize',14);
zlabel(names{p(3)},'FontSize',14);
grid
axis square
set(h,'PaperUnits','inches');
set(h,'PaperPosition',[.01 .01 3 3]);
% view(3);
% % % print('-dpng','-r600',['C:\Reza\ECGChannels',num2str(p(1)),num2str(p(2)),num2str(p(3)),'.png']);
% % % print('-deps','-r600',['C:\Reza\ECGChannels',num2str(p(1)),num2str(p(2)),num2str(p(3)),'.eps']);
