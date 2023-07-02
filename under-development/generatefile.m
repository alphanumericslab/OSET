clear;
close all;

addpath('E:\Programs\Baseline Wander Comparison\ECG Normal Sinus\30SecSegmentsParams');
addpath('E:\Programs\Baseline Wander Comparison\ECG Normal Sinus\30SecSegments');

for k = 0:190,
    filename = ['seg_' , num2str(k,'%3.3d')];
    load(filename,'data','chname');

    filename = ['seg_' , num2str(k,'%3.3d'),'_params'];
    load(filename);
    
    params{1}.mean = ECGmn;
    params{1}.std = ECGsd;

    params{1}.a = OptParams(1:end/3);
    params{1}.b = OptParams(end/3+1:2*end/3);
    params{1}.theta = OptParams(2*end/3+1:end);
    params{1}.label = chname(1:end-4);

    filename = ['C:\Documents and Settings\a\Desktop\ECG Parameter Database\Normal Sinus Database\','params_' , num2str(k,'%5.5d')];
    save(filename,'params');
end
