function [x fs mode] = LoadSeizureEEGFull2(path, fname)

% load data
load([path fname]);
reg = regexp(fname, '_');
varname = [fname(reg(2)+1:reg(4)-1) '_' num2str(str2double(fname(end-7:end-4)))];
md = fname(reg(2)+1:reg(3)-1);
if(isequal(md, 'interictal'))
    mode = 1;
elseif(isequal(md, 'preictal'))
    mode = 2;
elseif(isequal(md, 'test'))
    mode = 3;
end
eval(['data = ' varname '; clear ' varname ';']);
x = data.data;
fs = data.sampling_frequency;
clear data;

