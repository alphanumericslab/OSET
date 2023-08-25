function x = ECGModel(params, phase)
% ECGModel has been deprecated. Use ecg_gen_from_phase instead.
warning('ECGModel has been deprecated. Use ecg_gen_from_phase instead.');
x = ecg_gen_from_phase(params, phase);