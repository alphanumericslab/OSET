function R = Rotate3Dvec(tetax,tetay,tetaz);

L = max([length(tetax) , length(tetay) , length(tetaz)]);

Rx = zeros(3,3,L);
Ry = zeros(3,3,L);
Rz = zeros(3,3,L);

[1 0 0 ; 0 cos(tetax) sin(tetax) ; 0 -sin(tetax) cos(tetax)];
Ry = [cos(tetay) 0 -sin(tetay) ; 0 1 0 ; sin(tetay) 0 cos(tetay)];
Rz = [cos(tetaz) sin(tetaz) 0; -sin(tetaz) cos(tetaz) 0 ; 0 0 1];

R = Rx*Ry*Rz;
