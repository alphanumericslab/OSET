function [params model er] = ECGOptimizeModel(x, indexes, meanphase)
tetai = meanphase(indexes);
alphai = 1*x(indexes);
bi = .01*ones(size(alphai));

% options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100);
options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'Display','off');
InitParams = [alphai bi tetai];

% params = nlinfit(meanphase,x,@ECGModel,InitParams,options);
params = lsqnonlin(@(InitParams) ECGModelError(InitParams,x,meanphase,0),InitParams,InitParams-2,InitParams+2,options);

% Model0 = ECGModelError(InitParams,x,meanphase,1);
model = ECGModelError(params,x,meanphase,1);
er = 100*mean((x-model).^2)/mean(x.^2);
