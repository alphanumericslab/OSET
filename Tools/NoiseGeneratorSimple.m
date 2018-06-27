function noise =  NoiseGeneratorSimple(noisetype, NoisePower, fs, N)
%
% Simplified ECG noise generator for simulation purposes
switch(noisetype)
    case 'WHITE'
        noise = sqrt(NoisePower)*randn(1, N);
    case 'MA'
        load MA.mat MA
        ma = MA(:,3)';
        if(fs ~= 360)
            ma = resample(ma, fs, 360);
        end
        ma = ma - LPFilter(ma, 1.0/fs);
        ma = ma(1:N);
        ma = (ma - mean(ma))/std(ma);
        noise = sqrt(NoisePower)*ma;
    case 'EM'
        load EM.mat EM;
        em = EM(:,3)';
        if(fs ~= 360)
            em = resample(em, fs, 360);
        end
        em = em - LPFilter(em, 1.0/fs);
        em = em(1:N);
        em = (em - mean(em))/std(em);
        noise = sqrt(NoisePower)*em;
end
