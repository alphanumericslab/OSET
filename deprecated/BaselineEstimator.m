function baseline = BaselineEstimator(x_raw, baseline_removal_method, params)
% BaselineEstimator has been deprecated. Use baseline_filter instead.
warning('BaselineEstimator has been deprecated. Use baseline_filter instead.');
baseline = baseline_filter(x_raw, baseline_removal_method, params);