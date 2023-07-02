function y = TikhonovRegularization(x, diff_order_or_filter_coefs, lambda)
% TikhonovRegularization has been deprecated. Use tikhonov_regularization instead.
warning('TikhonovRegularization has been deprecated. Use tikhonov_regularization instead.');
y = tikhonov_regularization(x, diff_order_or_filter_coefs, lambda);