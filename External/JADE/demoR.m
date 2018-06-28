disp('This code is just a front-end to check that the real implementation');
disp('of JADE does run as seen on TV.');
disp(' ');
disp('This is indeed just a toy problem... Please run JADE on your real data.');
disp(' ');
disp(' ');
disp('The demo will display a figure with 3 rows: ');
disp('o First row: three `source signals'': a sine wave, a square wave, a');
disp('  whiste Gaussian process.');
disp('o Second row: the result of mixing by a 4x3 matrix and adding a');
disp('  small Gaussian noise.');
disp('o Thrid row: the result of Jade processing');
disp(' ');
disp(' ');
disp('The program repeatedly tries random mixing matrices.');
disp(' ');
disp(' ');
disp(' ');
disp('Hit <Return> key when ready.'); pause;
%
%
%
%Comments, bug reports, info requests are appreciated
%and should be directed to cardoso@sig.enst.fr (Jean-Francois Cardoso)
%Author : Jean-Francois Cardoso CNRS URA 820 / GdR TdSI / Telecom Paris

%=======================================================================
n	= 3 	;  % M = number of sources
m	= 4	;  % m = number of sensors (add the relevant lines in S= ...)
T	= 200	;  % sample size
NdB	= -30 	;  % kind of noise level in dB
%----------------------------------------------------------------------------

f1 = 0.013 ;
f2 = 0.02 ;

s1	=      cos(2*pi*f1*(1:T)) ;
s2	= sign(cos(2*pi*f2*(1:T))) ;
s3	=      randn(1,T) ;
s4	= sign(randn(1,T)) ;


S	= [ s1 ; s2 ; s3  ] ;


figure; clf 


disp('____________________________________________________________');


while 1


%%%%%%%%%%%  Source signals  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for is=1:n,
 	subplot(3,m,is);
	plot(S(is,:));
	axis([ 1 T -2 2 ]);
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('First row in figure: the source signals\n');

%%%%%%%%%%%  Mixing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mixing and noising


% random mixing matrix
A	= randn(m,n) ;
disp('Mixing matrix');
disp(A);


noiseamp 	= 10^(NdB/20) ;
X		= A*S + noiseamp*randn(m,T) ;
for ic=1:m
 	subplot(3,m,ic+m);
	plot(X(ic,:));
%	title('Hist. of observations');
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('Second row in figure: the observed mixtures\n');


%%%%%%%%% Separation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separation


fprintf('\nStrike any key to unmix\n');pause;
fprintf('\nIdentification running ......\n');
B 	= jadeR(X,n);
fprintf('\nIdentification completed ......\n');
Se	= B * X ;

for is=1:n,
 	subplot(3,m,is+m+m);
 	plot(Se(is,:));
%	title('Sep. histogram');
	axis([ 1 T -2 2 ]);
	set(gca,'Xtick',[1 T ]);
	set(gca,'Ytick',[]);
end;
drawnow; fprintf('Third row in figure: the estimated source signals\n');

% Performance
disp(' ');
disp(' ');
disp('Global system:');
disp('If this matrix is close to a product Permutation*Diagonal,')
disp('then separation was successful.');
% The so called `rejection rates'
disp((B*A).^2);

disp('____________________________________________________________');

fprintf('\nHit <Return> for another experiment with a different mixture\n');pause;
clf ;


end; %endless loop

