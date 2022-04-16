function [theta_min theta_max theta_mean theta_med] = ColumnSubspaceAngle(A, B)
% Reza Sameni
% December 2015

if(size(A, 1) ~= size(B, 1))
   error('A and B should have the same number of rows'); 
end

ang = zeros(size(A, 2), size(B, 2));
for i = 1 : size(A, 2),
    for j = 1 : size(B, 2),
        ang(i, j) = 90*acos(A(:,i)'*B(:,j)/sqrt((A(:,i)'*A(:,i)) * (B(:,j)'*B(:,j))))/(pi/2);
    end
end

theta_min = min(ang(:));
theta_max = max(ang(:));
theta_mean = mean(ang(:));
theta_med = median(ang(:));