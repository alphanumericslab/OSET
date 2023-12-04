% This code is just a front-end to source separation algorithms.
% Purpose:
% 1) generate synthetic data
% 2) call some source separation algorithm
% 3) display the results
% The data are CM (constant modulus signals and QAM4.
% The mixing matrix  is randomly generated.
% Comments, bug reports, info requests are appreciated
% and should be directed to cardoso@sig.enst.fr (Jean-Francois Cardoso)
% Author : Jean-Francois Cardoso CNRS URA 820 / GdR TdSI / Telecom Paris

%=======================================================================
N	= 4	;  % N = number of sensors (add the relevant lines in S= ...)
M	= 3 	;  % M = number of sources
T	= 200	;  % sample size
NdB	= -15 	;  % kind of noise level in dB
%----------------------------------------------------------------------------
disp('Each of the following plots shows the COMPLEX PLANE.')
disp('Each point is a sample of a source signal, of a sensor output')
disp('or a separated signal as indicated.')


while 1
disp('____________________________________________________________');

% the source signals
S= [ ...
exp(2*i*pi*rand(1,T))					; 	% constant modulus random phase
exp(2*i*pi*rand(1,T))					; 	% constant modulus random phase
(2*fix(2*rand(1,T))-1+i*(2*fix(2*rand(1,T))-1))/sqrt(2) ;	% QAM4
];

% random mixing matrix
A=randn(N,M)+j*randn(N,M);
disp('Mixing matrix');disp(A);

clf;
subplot(1,1,1);
for is=1:M,
 subplot(2,2,is);
 plot(S(is,:),'.');title('One of the source signals');
 axis('square'); axis('equal'); axis([-2 2 -2 2]);
end;
fprintf('\nStrike any key to mix\n');pause;


% mixing and noising
noiseamp = 10^(NdB/20)/sqrt(2) ; % (the sqrt(2) accounts for real+imaginary powers)
X= A*S + noiseamp*(randn(N,T)+i*randn(N,T));

clf;
for is=1:min([ N 4]),
 subplot(2,2,is);
 plot(X(is,:),'.');title('One of the mixed signals');
 axis('square');axis('equal');%axis([-2 2 -2 2]);
end;
fprintf('\nStrike any key to unmix\n');pause;

% Separation
fprintf('\nIdentification running ......\n');
[Ae,Se]=jade(X,M);

clf;
for is=1:M,
 subplot(2,2,is);
 plot(Se(is,:),'.');title('One of the separated signals');
 axis('square');axis('equal');axis([-2 2 -2 2]);
end;

% Performance
disp('Performance:');
disp('The global (sepration*mixing) matrix should be close to a permutation');
disp('The following shows its squared entries (rejection levels)');
disp(abs(pinv(Ae)*A).^2);

fprintf('\nStrike any key for a different mixture\n');pause;

end; %endless loop

