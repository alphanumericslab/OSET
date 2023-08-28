function [dipole, teta] = DipoleGenerator3(HR, fs, alphai, bi, tetai, teta0)
% DipoleGenerator3 has been deprecated. Use vcg_gen_var_hr instead with additional input parameter method = 'LOOP'.
warning('DipoleGenerator3 has been deprecated. Use vcg_gen_var_hr instead with additional input parameter method = ''LOOP''');
method = 'LOOP';
[dipole, teta] = vcg_gen_var_hr(HR, fs, alphai, bi, tetai, teta0, method);