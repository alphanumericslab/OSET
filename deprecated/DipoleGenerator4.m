function [dipole, teta] = DipoleGenerator4(HR, fs, alphai, bi, tetai, teta0)
% DipoleGenerator4 has been deprecated. Use dipole_gen_var_hr instead with additional input parameter method = 'MATRIX'.
warning('DipoleGenerator4 has been deprecated. Use dipole_gen_var_hr instead with additional input parameter method = ''MATRIX''');
method = 'MATRIX';
[dipole, teta] = dipole_gen_var_hr(HR, fs, alphai, bi, tetai, teta0, method);