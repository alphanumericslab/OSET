function E = ECGModelError(X,ECGmn,Phasemn,flag)
% ECGModelError has been deprecated. Use ecg_gen_from_phase instead.
warning('ECGModelError has been deprecated. Use ecg_gen_from_phase instead with required changes.');

Z = ecg_gen_from_phase(X, Phasemn);
if flag==0
    E = Z - ECGmn;
elseif flag==1
    E = Z;
end