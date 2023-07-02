% test implementation of periodic splines
% UNDER TEST
% Reza Sameni, Copyright 2015
%

close all;
clear;
clc;

fsim = 1000; % Simulation sampling frequency (Hz)
knots = [1 3 5 7 15 25 35]; % knot points (s)
D = 4; % degree
tfirst = 0; % s
tlast = 25; % s

N = length(knots) - 1;
if(D > (N-1))
    warning('D may not exceed (N-1). Reduce D, or add some knots!');
end

delta = 1/fsim;

A = zeros( N*(D + 1) - D*(D + 1)/2 ); % calculated to be this size!
row = 0;
% first_row_in_D = zeros(1, D + 1);
B0 = zeros(size(A, 1), 1);
for k = 0 : D,
    for i = 0 : N - k - 1,
        row = row + 1;
        if(k > 0)
            col1 = first_row_in_previous_level + i;
            col2 = first_row_in_previous_level + i + 1;
            A(row, col1) = k/(knots(mod(i + k, N + 1) + 1) - knots(mod(i, N + 1) + 1));
            A(row, col2) = -k/(knots(mod(i + k + 1, N + 1) + 1) - knots(mod(i + 1, N + 1) + 1));
            B0(row) = (tfirst - knots(mod(i, N + 1) + 1))/(knots(mod(i + k, N + 1) + 1) - knots(mod(i, N + 1) + 1)) + ...
                      (knots(mod(i + k + 1, N + 1) + 1) - tfirst)/(knots(mod(i + k + 1, N + 1) + 1) - knots(mod(i + 1, N + 1) + 1));
        else
            if(tfirst >= knots(mod(i, N + 1) + 1) && tfirst < knots(mod(i + 1, N + 1) + 1))
                B0(row) = 1;
            else
                B0(row) = 0;
            end
        end
        if(i == N - k - 1)
            first_row_in_previous_level = row - (N - k - 1);
        end
    end
end

AA = delta * A + eye(size(A));

% t = linspace(min(knots), max(knots), round((max(knots) - min(knots))*fsim));
t = linspace(tfirst, tlast, round((tlast - tfirst)*fsim));

B = zeros(size(A,1), length(t));
B(:, 1) = B0;
for i = 2:size(B, 2),
    B(:, i) = AA * B(:, i - 1);
end

figure
plot(t, B');
grid
