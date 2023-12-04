clear;
close all;

% Synthetic maternal parameters
params{1}.mean = [];
params{1}.std = [];

params{2}.mean = [];
params{2}.std = [];

params{3}.mean = [];
params{3}.std = [];

params{1}.theta  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68];
params{1}.a = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39];
params{1}.b     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020    0.1812];

params{2}.theta  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8      1.58];
params{2}.a = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08];
params{2}.b     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3];

params{3}.theta  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55];
params{3}.a = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35];
params{3}.b     = [.03  .12  .04         .4    .045       .05    .8 .4 .2];

params{1}.label = 'x-axis of synthetic maternal VCG';
params{2}.label = 'y-axis of synthetic maternal VCG';
params{3}.label = 'z-axis of synthetic maternal VCG';

save('C:\Documents and Settings\a\Desktop\ECG Parameter Database\Misc\params_02000.mat','params');

%//////////////////////////////////////////////////////////////////////////
% % Synthetic fetal parameters
params{1}.mean = [];
params{1}.std = [];

params{2}.mean = [];
params{2}.std = [];

params{3}.mean = [];
params{3}.std = [];

params{1}.theta  = [-0.7    -0.17    0       0.18     1.4];
params{1}.a = [0.07     -0.11   1.3     0.07   0.275];
params{1}.b     = [.1       .03     .045     0.02    0.3];

params{2}.theta  = [-0.9     -0.08   0       0.05        1.3];
params{2}.a = [0.04     0.3     .45     -0.35       0.05];
params{2}.b     = [.1       .05      .03    .04         .3];

params{3}.theta  = [-0.8      -.3     -0.1        .06     1.35];
params{3}.a = [-0.14    .03     -0.4        .46     -0.1];
params{3}.b     = [.1       .4      .03         .03     .3];

params{1}.label = 'x-axis of synthetic fetal VCG';
params{2}.label = 'y-axis of synthetic fetal VCG';
params{3}.label = 'z-axis of synthetic fetal VCG';

save('C:\Documents and Settings\a\Desktop\ECG Parameter Database\Misc\params_02001.mat','params');