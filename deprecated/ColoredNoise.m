function n = ColoredNoise(sd,Len,fs,beta)
% ColoredNoise has been deprecated. Use colored_noise_gen instead.
warning('ColoredNoise has been deprecated. Use colored_noise_gen instead.');
n = colored_noise_gen(sd,Len,fs,beta);