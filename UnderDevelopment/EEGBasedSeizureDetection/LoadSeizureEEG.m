function [x fs sequence] = LoadSeizureEEG(path, subject, mode, number, ch, time)

load([path subject '_' mode '_segment_' num2str(number,'%.4d') '.mat']);
varname = [mode '_segment_' num2str(number,'%d')];
eval(['data = ' varname '; clear ' varname ';']);

x = data.data;
fs = data.sampling_frequency;
sequence = data.sequence;

samples = round(fs*time);
samples(samples < 1) = 1;
samples(samples > size(x,2)) = size(x,2);
x = x(ch, samples(1):samples(2));

