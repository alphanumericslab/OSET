function [A,B,ang] = RandomMatrices(angle,dim1,flag)
% [A,B,ang] = RandomMatrices(angle,dim,flag),
% 
% A random matrix generator with a specific degree between its column
% space.
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.



Kmax = 10000;

k = 0;
if flag == 0,
    ang = [90 90 90];
    while max(ang)>angle && k < Kmax,
        A = randn(dim1,3);
        R = [diag(rand(dim1-3,1))*eye(dim1-3,dim1) ; zeros(3,dim1-3) Rotation3D(pi*angle*rand/180,pi*angle*rand/180,pi*angle*rand/180)];

        B = R*A;

        [Qa,Ra] = qr(A,0);
        [Qb,Rb] = qr(B,0);
        C = Qa'*Qb;
        [U,S,V] = svd(C);
        S = min(diag(S),1);
        S = max(S,-1);
        ang = acos(S)*180/pi;

        k = k+1;
    end

elseif flag ==1,
    ang = [0 0 0];
    while min(ang)<angle && k < Kmax,
        A = randn(dim1,3);
        B = randn(dim1,3);

        [Qa,Ra] = qr(A,0);
        [Qb,Rb] = qr(B,0);
        C = Qa'*Qb;
        [U,S,V] = svd(C);
        S = min(diag(S),1);
        S = max(S,-1);
        ang = acos(S)*180/pi;

        k = k+1;
    end

end

if(k>=Kmax)
    A = [];
    B = [];
    ang = [];
end

%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////
function R = Rotation3D(tetax,tetay,tetaz)
Rx = [1 0 0 ; 0 cos(tetax) sin(tetax) ; 0 -sin(tetax) cos(tetax)];
Ry = [cos(tetay) 0 -sin(tetay) ; 0 1 0 ; sin(tetay) 0 cos(tetay)];
Rz = [cos(tetaz) sin(tetaz) 0; -sin(tetaz) cos(tetaz) 0 ; 0 0 1];
R = Rx*Ry*Rz;
