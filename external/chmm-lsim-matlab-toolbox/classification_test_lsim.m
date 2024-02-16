addpath(genpath(pwd))
clear
close all
clc


mydir = pwd;
idcs = strfind(mydir,filesep);
% second parent folder contains the datasets
dataset_dir = [mydir(1:idcs(end-1)-1),'/DataSets'];
results_dir = [mydir(1:idcs(end-1)-1),'/Results/',mydir(idcs(end-1)+1:end)];
results_dir = pwd;
addpath([mydir(1:idcs(end-1)-1),'/chmm-lsim-karimi-toolbox'])

%%

test_number = 10;
max_repeat = 1000;


for C = 3 %C is number of channels in CHMM

    for T = 5:5:10    % T  number of time samples

        k1 = 0;
        for train_number = [10,15,20]

            k1=k1+1;

            parfor repeat_num = 1:max_repeat
                clc
                C
                T
                train_number
                [scores_lsim(repeat_num,:), Label(repeat_num,:)] = sim_calassification_lsim(C,T,train_number, test_number);
            end

            %             Label = repmat(Label,max_repeat,1);
            ACC_lsim(C,T,train_number) = sum((scores_lsim(:)>0 &Label(:)>0.5)|(scores_lsim(:)<=0 & Label(:)<0.5))/(length(Label(:)));

            [~,~,~,AUC] = perfcurve(Label(:),scores_lsim(:),+1);
            AUC_lsim(C,T,train_number) = AUC;
            close all

            save([results_dir,'/lsim.mat'])

        end

    end
end


save([results_dir,'/lsim.mat'])


%% Plot

clc
load([results_dir,'/lsim.mat'])


C = 2; %C is number of channels in CHMM
T = 5:5:10;    % T  number of time samples
train_number = [10,15,20];

ACC_lsim(C,T,train_number) 




