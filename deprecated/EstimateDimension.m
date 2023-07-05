function [lambda, AIC, MDL, NEW, ENSTh, ENS, fhand] = EstimateDimension(x, noisevar, flagplot)
% EstimateDimension has been deprecated. Use dimension_estimation instead.
warning('EstimateDimension has been deprecated. Use dimension_estimation instead.');
[lambda, AIC, MDL, NEW, ENSTh, ENS, fhand] = dimension_estimation(x, noisevar, flagplot);