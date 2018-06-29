% A test of correlation dimension calculation
%

clear;
close all;

fs = 1000;
data = load('s0010_re.txt');
% data = load('s0014lre.txt');


time = data(:,1)';
data = data(:,2:end)';

data = data - LPFilter(data,1.5/fs);
data = LPFilter(data,150/fs);

I = 1:2*fs;
data = data(:,I);

% data = randn(size(data,1),5000);
data = .41*randn(size(data));


N = size(data,2);

r = .0001:.1:2*max(max(abs(data)));
% r = log(.0001):.1:log(max(max(abs(data))));
% r = exp(r);

C = zeros(1,length(r));
for k = 1:length(r)
    c = 0;
    for i = 1:N-1
        Xj = data(:,i+1:N);
        Xi = data(:,i(ones(size(Xj,2),1)));
        d = sqrt(sum((Xi-Xj).^2,1));
        n = find(d<r(k));
        c = c + length(n);
    end
    C(k) = c*2/N/(N-1);
end

figure
plot(log(r),log(C),'ro');
grid;
% axis tight
xlabel('log(r)','FontSize',16);
ylabel('log(C)','FontSize',16);
set(gca,'Box','On','FontSize',16);
set(gcf,'Position',[520 678 580 420]);

% % % figure
% % % % plot(r,C);
% % % % hold on;
% % % plot(r,C,'ro');
% % % grid;
% % % xlabel('r');
% % % ylabel('C');

slope = diff(log(C))./diff(log(r));
figure;
% plot(r(1:end-1),slope);
% hold on;
plot(r(1:end-1),slope,'ro','linewidth',2);
grid
axis tight
xlabel('r','FontSize',16);
ylabel('$\frac{\partial log C(r)}{\partial log(r)}$','FontSize',16,'Interpreter','latex');
set(gca,'Box','On','FontSize',16);

% % % plot(log(r),log(C));
% % % hold on;
% % % plot(log(r),log(C),'r.');
% % % grid;

% % % % % 
% % % % % 
% % % % % W = jadeR(data);
% % % % % s = W*data;
% % % % % 
% % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % I = 1:10*fs;
% % % % % L = 5;
% % % % % k = 0;
% % % % % names = {'i','ii','iii','avr','avl','avf','v_1','v_2','v_3','v_4','v_5','v_6','v_x','v_y','v_z'};
% % % % % for i = 1:size(s,1),
% % % % %     if(mod(i,L)==1)
% % % % %         h = figure;
% % % % %         k = k + 1;
% % % % %     end
% % % % %     subplot(L,1,mod(i-1,L)+1);
% % % % %     plot(time(I),data(i,I),'k');
% % % % %     
% % % % %     a = axis;
% % % % %     a(1) = time(1);
% % % % %     a(2) = 10;
% % % % %     axis(a);
% % % % %     ylabel(names{i},'FontSize',16);
% % % % %     grid;
% % % % %     set(gca,'Box','On','FontSize',16);
% % % % %     if (mod(i,L)~=0)
% % % % %         set(gca,'XTickLabel',[]);
% % % % %     else
% % % % %         xlabel('Time (s)','FontSize',16);
% % % % %     end
% % % % % 
% % % % %     if(mod(i,L)==0)
% % % % %         set(h,'PaperUnits','inches');
% % % % %         set(h,'PaperPosition',[.01 .01 2.5 8.5])
% % % % % % % %         print('-dpng','-r600',['C:\Reza\ECGChannels_',num2str(k),'.png']);
% % % % % % % %         print('-deps','-r600',['C:\Reza\ECGChannels_',num2str(k),'.eps']);
% % % % %     end
% % % % % end
% % % % % 
% % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % I = 1:10*fs;
% % % % % L = 5;
% % % % % k = 0;
% % % % % for i = 1:size(s,1),
% % % % %     if(mod(i,L)==1)
% % % % %         h = figure;
% % % % %         k = k + 1;
% % % % %     end
% % % % %     subplot(L,1,mod(i-1,L)+1);
% % % % %     plot(time(I),s(i,I),'k');
% % % % %     
% % % % %     a = axis;
% % % % %     a(1) = time(1);
% % % % %     a(2) = 10;
% % % % %     axis(a);
% % % % %     ylabel(['IC_{',num2str(i),'}'],'FontSize',16);
% % % % %     grid;
% % % % %     set(gca,'Box','On','FontSize',16);
% % % % %     if (mod(i,L)~=0)
% % % % %         set(gca,'XTickLabel',[]);
% % % % %     else
% % % % %         xlabel('Time (s)','FontSize',16);
% % % % %     end
% % % % % 
% % % % %     if(mod(i,L)==0)
% % % % %         set(h,'PaperUnits','inches');
% % % % %         set(h,'PaperPosition',[.01 .01 2.5 8.5])
% % % % % % % %         print('-dpng','-r600',['C:\Reza\ICChannels_',num2str(k),'.png']);
% % % % % % % %         print('-deps','-r600',['C:\Reza\ICChannels_',num2str(k),'.eps']);
% % % % %     end
% % % % % end
% % % % % %//////////////////////////////////////////////////////////////////////////
% % % % % 
% % % % % p = [7 8 11];
% % % % % h = figure;
% % % % % plot3(data(p(1),I),data(p(2),I),data(p(3),I),'b.');
% % % % % set(gca,'Box','On','FontSize',14);
% % % % % xlabel(names{p(1)},'FontSize',14);
% % % % % ylabel(names{p(2)},'FontSize',14);
% % % % % zlabel(names{p(3)},'FontSize',14);
% % % % % grid
% % % % % axis square
% % % % % set(h,'PaperUnits','inches');
% % % % % set(h,'PaperPosition',[.01 .01 3 3]);
% % % % % % view(3);
% % % % % print('-dpng','-r600',['C:\Reza\ECGChannels',num2str(p(1)),num2str(p(2)),num2str(p(3)),'.png']);
% % % % % print('-deps','-r600',['C:\Reza\ECGChannels',num2str(p(1)),num2str(p(2)),num2str(p(3)),'.eps']);


% % % L1 = 10;
% % % figure;
% % % for i = 1:L1,
% % %     subplot(L1,1,i);
% % %     plot(time,s(i,:),'b','linewidth',1);
% % %     grid;
% % %     set(gca,'Box','On','FontSize',16);
% % %     if (i<L1)
% % %         set(gca,'XTickLabel',[]);
% % %     end
% % %     ylabel(['IC_',num2str(i)],'FontSize',16);
% % %     if (i==1),
% % %         title('Fetal MCG After Maternal MCG Removal and Decomposed into Periodic Components (using \PiCA)');
% % %     end
% % % end


% PlotECG(s,4,'b');