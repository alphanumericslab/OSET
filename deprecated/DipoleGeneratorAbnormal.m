function [dipole, teta]= DipoleGeneratorAbnormal(N,fs,f,alphai,bi,tetai,teta0,STM,S0)
% DipoleGeneratorAbnormal has been deprecated. Use vcg_gen_abnormal instead.
warning('DipoleGeneratorAbnormal has been deprecated. Use vcg_gen_abnormal instead.');
[dipole, teta]= vcg_gen_abnormal(N,fs,f,alphai,bi,tetai,teta0,STM,S0);