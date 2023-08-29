new_path = '/Users/rsameni/Documents/GitHub/OSET/external/MutualInformationICA';
current_path = getenv('PATH');
setenv('PATH', [new_path ':' current_path]);

%The first step is to compile the C-source codes into executable files, which will be then called be the MATLAB Programms
% for example  (under linux)
%    icc -c miutils.C -o miutils.o 
%    icc  miica.C -o milca miutils.o 

% When you have now the executable file (e.g. milca) you can call already the algorithm.
%   for example   
%   [icasignal, W_matrix]=milca(input_signal);


% Here are now examples for each Programm

% MI Values -----------------------------------------------------------------
fprintf('\n');
fprintf('Examples: MI Calculations in different dimensions and spaces \n');
fprintf('\n');
fprintf('2 independent uniformly distributed signals x \n');
x=rand(2,3000);

fprintf('MI should be around 0 (because of statistical fluctuations it can produce also small negativ values)\n');
MIhigherdim(x)

fprintf('mix the two signals - make them dependent:  y=A*x\n');
A=rand(2);
y=A*x;

fprintf('MI of the mixed signal:  \n');
MIhigherdim(y)

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('the same is valid for higher dimensions\n');

x=rand(6,3000);
fprintf('MI of the 6-dimensional signal:\n');
MIhigherdim(x)
A=rand(6);
y=A*x;
fprintf('MI of the mixed 6-dimensional signal:\n');
MIhigherdim(y,6,1,1)

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('we embedded the signal now to inlcude also dependencies in time (time structure)\n');
fprintf('Because in our case there is no time structure in the signal the obtained MI should be similar (because of the higher dimension a little bit higher)\n');
MIhigherdim(y,6,2,2)

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('Now we have a 8 dimensional signal where only the first 4 and last 4 components are dependent among each other\n');

x=rand(8,2000);
A=rand(4);
B=rand(4);
C=eye(8);
C(1:4,1:4)=A;
C(5:8,5:8)=B;
y=C*x;


fprintf('The MI of the 8 dimensional signal should be high: \n');
MIhigherdim(y)
fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;
fprintf('But the MI between x=channel 1-4  and y= 5-8 should be again around zero\n');
MIxnyn(y(1:4,:),y(5:8,:))
fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;
fprintf('The MI between x=channel 1-3  and y= 4-8 should be higher than zero because channel 4 contains information of x\n');
MIxnyn(y(1:3,:),y(4:8,:))

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

%MILCA --------------------------------------------------------------------------


fprintf('\n');
fprintf('Using the accurate MI estimator as a contrast function in the ICA algorithm: \n');
fprintf('\n');
fprintf('Generate 3 dimensional signal with 1 uniformly distributed signal,\n');
fprintf('1 white Gaussian, and 1 red Gaussian');
[B,A] = butter(6,0.3);
x=[rand(1,2000);randn(1,2000);filter(B,A,randn(1,2000))];
A=rand(3);
fprintf('mix the three signals:  mix=A*x\n');
mix=A*x;

fprintf('Demix signal with milca  -  [icasignal, W]=milca(mix) ');
[icasignal, W]=milca(mix);    

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('Reliability of the ICA output in the first figure (diagonale set to zero) -\n')
fprintf('dependency matrix (MI between all channel combination) one can see that all ICA components\n');
fprintf('are independent from each other - the variability matrix (diagonale set to zero) show \n');
fprintf('that 2 components are not reliable separable - low value (because the both have Gaussian distribution), but ... \n');

ICAtests(icasignal,1);

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('because of the the different spectra of the Gaussian signals one can use\n');
fprintf('also the time information contained in the signal\n');
fprintf('The same tests with embedded signals show that indeed the two Gaussian\n');
fprintf('signal are dependent (figure 1 (left) and figure 2 - minimum not by angle zero)\n'); 
% P.S. with a very low probability it can hapen that the correct de-mixing 
% matrix (angle) is estimated by chance (simply run algorithm again)

ICAtests(icasignal,2);

fprintf('\n');
fprintf('Pause - press any key to continue\n');
fprintf('\n');
pause;

fprintf('Now we use also the time information in the ICA-algorithm - [icasignaldelay, W]=milcadelay(mix)\n');
fprintf('please wait (1-2 minutes)\n');
fprintf('\n');
[icasignaldelay, W]=milcadelay(mix); 

fprintf('The tests shows that all components are now independent (also in time)\n');
fprintf('and all components are reliable estimates (variability matrix - only high values)\n'); 

ICAtests(icasignaldelay,2);

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

%Clustering --------------------------------------------------------------------------


fprintf('\n');
fprintf('Using the MI estimator to cluster dependent data: \n');
fprintf('\n');
fprintf('Load an ICA output of a 8 channel ECG of a pregnant woman\n');
% (www.esat.kuleuven.ac.be/sista/daisy)
load milcaECG8ch

fprintf('The ICA components are still dependent \n');
ICAtests(icasig,1);

fprintf('\n');
fprintf('Pause - press any key to continue\n');
pause;

fprintf('we apply hierarchical clustering using the grouping property of mutual\n');
fprintf('information to obtain a dendrogram of the dependencies.\n');
fprintf('Uses the MATLAB function dendrogram, which is included only in the statistical toolbox\n');
MIClustering(icasig);

fprintf('If you want to improve this result look into the Matlab file (Tutorial.m) \n');


% % Because we are not completely satisfied with the result we proceed like
% % in [1] - embed the 8 channels in 24 dimensions and apply then MILCA
% [Nd,N]=size(ecgsignal);
% tau=2;
% y1=ecgsignal(:,1:(N-2*tau));
% y2=ecgsignal(:,(tau+1):(N-tau));
% y3=ecgsignal(:,(2*tau+1):N);
% ecgsignal=[y1;y2;y3];
% 
% This can take same time, depending on our computer (~ hour)
% [ica24, W]=milca(ecgsignal); 
% 
% % We can separate now the mother and child contributions better
% ICAtests(ica24,1);
% MIClustering(ica24);
%
% % to see the difference between 8 and 24 channels make a back-projection
% % only of the mother, foetus cluster, respectively (see [1] Fig 20).


