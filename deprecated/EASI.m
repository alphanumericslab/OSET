function [y, B, conv] = EASI(x, nsou, lambda, nlintype)
% EASI has been deprecated. Use easi_source_separation instead.
warning('EASI has been deprecated. Use easi_source_separation instead.');
[y, B, conv] = easi_source_separation(x, nsou, lambda, nlintype);