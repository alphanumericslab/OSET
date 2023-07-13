function [ARCoefsKF , ARCoefsKS] = TimeVariantAR(data, order, x0, q, R, p0, alpha)
% TimeVariantAR has been deprecated. Use time_variant_ar_kalman instead.
warning('TimeVariantAR has been deprecated. Use time_variant_ar_kalman instead.');
[ARCoefsKF , ARCoefsKS] = time_variant_ar_kalman(data, order, x0, q, R, p0, alpha);