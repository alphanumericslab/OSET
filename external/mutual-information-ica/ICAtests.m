function ICAtests(x,emb_dim,emb_tau,kneig,algo,kmax)

% Reliability test for any ICA output (linear transformation version!)
% Output: Dependency Matrix, Variability Matrix, and MI vs rotation angle
% plot (to check if the components are indeed the most independent one)

% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm
% algo ... 1=cubic  2=rectangular
% kmax ... Number of angles between 0 and pi/2 
% emb_dim... embedding dimension (default is 1, no embedding)
% emb_tau... time-delay, only relevant when emb_dim>1 (default is 1)

%default-values
if ~exist('emb_dim'), emb_dim=1; end
if ~exist('emb_tau'), emb_tau=1; end
if ~exist('kneig'), kneig=5; end
if ~exist('algo'), algo=2; end
if ~exist('kmax'), kmax=15; end


[Nd,N]=size(x);
if Nd>N
    x=x';
    [Nd,N]=size(x);
end

% save data for external Programm
zwsp=x';
save zwsptests.txt zwsp -ASCII

% execute C Programm
[a b]=unix(['ICAtests zwsptests.txt ',num2str(Nd),' ',num2str(N),' ',num2str(kneig),' ',num2str(emb_tau),' ',num2str(emb_dim),' ',num2str(algo-1),' ',num2str(kmax)]);

%format output
[mivsangle,count,errmsg,nextindex] = sscanf(b,'%f',[kmax,Nd*(Nd-1)/2]);
varmimat=str2num(b(nextindex:end));
mimat=varmimat(1:Nd,:);
varmimat=varmimat(Nd+1:2*Nd,:);

  
%plot output
count=1;
mivsangle=mivsangle';

if Nd==2   
    fprintf('MI value: %f \n',mimat(1,2));
    fprintf('Variability of the MI under rotation: %f',varmimat(1,2));
    figure;
    plot(mivsangle);
    set(gca,'xlim',[1 kmax],'xtick',round(linspace(1,kmax,4)),'xticklabel',(round(linspace(1,kmax,4))-1)/(kmax)*90);
    xlabel('roation angle');
    ylabel('mutual information');
else
    figure;
    subplot(1,2,1)
    maxMI=max(max(mimat));
    if maxMI<0.3, maxMI=0.3; end  
    imagesc(mimat,[0,maxMI]);
    title('Dependency Matrix');
    h=flipud(colormap(gray));
    colormap(h);
    colorbar;
    subplot(1,2,2);
    minVar=min(min(varmimat));
    if minVar<-0.015, fprintf('Input seems to be not the most independent representation under linear transformation !!!!!!!!!!!'); end
    imagesc(varmimat);
    title('Variability Matrix');
    colormap(h);
    colorbar;
    
    figure;
    
    ymin=min(min(mivsangle));
    ymax=max(max(mivsangle));
    for i=1:Nd
        for j=1:Nd
            if i<j
                subplot(Nd,Nd,j+(i-1)*Nd);
                plot(mivsangle(count,:),'.-');
                set(gca,'xlim',[1 kmax],'ylim',[ymin ymax],'xtick',round(linspace(1,kmax,4)),'xticklabel',(round(linspace(1,kmax,4))-1)/(kmax)*90);
                if count==1
                    xlabel('roation angle');
                    ylabel('mutual information');
                end
                count=count+1;
            end
        end
    end
end





