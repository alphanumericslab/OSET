function [x fs] = LoadSeizureEEGFull(path, subject, mode, number)

load([path subject '_' mode '_segment_' num2str(number,'%.4d') '.mat']);
varname = [mode '_segment_' num2str(number,'%d')];
eval(['data = ' varname '; clear ' varname ';']);

x = data.data;
fs = data.sampling_frequency;

