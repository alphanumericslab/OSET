% simple preprocessing of the raw data using OSET tools
clear;
close all;

for i = 1:20,
    inputfile = ['signal_',num2str(i,'%2.2d'),'.mat'];
    outputfile = ['signal_',num2str(i,'%2.2d'),'_preprocessed.mat'];
    load(inputfile);
    [num,den] = iirnotch(50/(fs/2),1/(fs/2));
    for j = 1:size(s,1),
        s(j,:) = filter(num,den,s(j,:));
%       b = TrimmedFilter(s(j,:),'trmean',(.4*fs),(.1*fs),(.1*fs));
%       b = TrimmedFilter(s(j,:),'median',(.4*fs),.01*fs,.2*fs);
        b = TrimmedFilter(s(j,:),'median',(.2*fs));
        b = TrimmedFilter(b,'median',(.4*fs));
        b = LPFilter(b,10/fs);
        s(j,:) = s(j,:) - b;
    end
    I = 1:20000;
    PlotECG(s(:,I),8,'b',fs);
    s = single(s);
    fs = single(fs);
    save(outputfile,'s','fs','-V6');
end