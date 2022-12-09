function baseline = BaselineEstimator(x_raw, baseline_removal_method, params)
% Different approaches for baseline wander estimation using OSET functions
%
% Reza Sameni, Nov 2021
%

switch baseline_removal_method
    case 'BYPASS' % bypass baseline wander (zero baseline)
        baseline = zeros(size(x_raw));
    case 'LP'
        baseline = LPFilter(x_raw, params.fc/params.fs);
    case 'MNMN' % two-stage moving average
        baseline = BaseLine2(x_raw, params.wlen1, params.wlen2, 'mn');
    case 'MDMN' % two-stage moving median plus moving average
        % apply a two-stage moving-median (md) and moving-average (mn) baseline detector
        bl1 = BaseLine1(x_raw, params.wlen1, 'md');
        baseline = BaseLine1(bl1, params.wlen2, 'mn');
    case 'MDMDMDMN' % three median filters followed by a moving average
        bl1 = BaseLine1(x_raw, params.wlen1, 'md');
        bl2 = BaseLine1(x_raw, params.wlen2, 'md');
        bl3 = BaseLine1(x_raw, params.wlen3, 'md');
        baseline = BaseLine1((bl1 + bl2 + bl3) / 3, params.wlen4, 'mn');
    case 'TIKHONOV' % Tikhonov regularization
        baseline = TikhonovRegularization(x_raw, params.DiffOrder, params.lambda);
    otherwise
        error('Unknown baseline removal method');
end