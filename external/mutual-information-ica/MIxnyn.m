function miout=MIxnyn(x,y,kneig);

% Calculate MI value between 2 vector of any dimension (rectangular
% version)
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm


%default-values
if ~exist('kneig'), kneig=6; end


% check input data if format is correct
[Ndx,Nx]=size(x);
if Ndx>Nx
    x=x';
    [Ndx,Nx]=size(x);
end
[Ndy,Ny]=size(y);
if Ndy>Ny
    y=y';
    [Ndy,Ny]=size(y);
end

if Nx~=Ny
    if Nx>Ny
        N=Ny;
    else
        N=Nx;
    end
    fprintf('The two input vectors must have the same length !!!!');
    fprintf('Caluculation using the %d datapoints',N);
    
else
    N=Nx;    
end


% save data for C-Programm
zwsp=[x;y]';
save zwspMIxnyn.txt zwsp -ASCII


% execute C Programm
[a unout]=unix(['MIxnyn zwspMIxnyn.txt ',num2str(Ndx),' ',num2str(Ndy),' ',num2str(N),' ',num2str(kneig)]);
miout=str2num(unout);


