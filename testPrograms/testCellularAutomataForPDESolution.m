% A simulation of the diffusion PDE numeric solution
% A comparison of two discretization methods
%
% Reza Sameni
% Apr 2021

clear
close all
clc

Lx = 51; % grid dimension 1
Ly = 51; % grid dimension 2
dx = 0.1; % in meters
dy = 0.1; % in meters

dt = 0.01; % simulation time unit
NT = 100; % total simulation time iterations

D = 1e-1; % diffusion parameter

C1 = zeros(Lx, Ly, NT); % method 1 concentration
C2 = zeros(Lx, Ly, NT); % method 2 concentration

% Ininitialization
C1(round(Lx/2), round(Ly/2), 1) = 1.0; % Point 1
C1(round(Lx/3), round(Ly/4), 1) = 0.7; % Point 2

C2(:, :, 1) = C1(:, :, 1); % both methods start with the same initial conditions

alpha_x = D * dt / dx^2;
alpha_y = D * dt / dy^2;
alpha_xy = D * dt / (dx^2 + dy^2);

% Check validity of the parameters
if((1 - 2 * alpha_x - 2 * alpha_y) < 0 || (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) < 0)
    error('Stability condition for parameters not fulfilled. Make simulation time period smaller');
end
    
h = figure;
for t = 1 : NT - 1
    % Update all points using two methods
    for i = 2 : Lx - 1
        for j = 2 : Ly - 1
            % 4-point neighborhood update
            C1(i, j, t + 1) = (1 - 2 * alpha_x - 2 * alpha_y) * C1(i, j, t) + alpha_x * C1(i - 1, j, t) + alpha_x * C1(i + 1, j, t) + alpha_y * C1(i, j - 1, t) + alpha_y * C1(i, j + 1, t);

            % 8-point neighborhood update
            C2(i, j, t + 1) = (1 - 2 * alpha_x - 2 * alpha_y - 2 * alpha_xy - 2 * alpha_xy) * C2(i, j, t) + alpha_x * C2(i - 1, j, t) + alpha_x * C2(i + 1, j, t) + alpha_y * C2(i, j - 1, t) + alpha_y * C2(i, j + 1, t) +   ...                 
                                alpha_xy * C2(i - 1, j - 1, t) + alpha_xy * C2(i + 1, j + 1, t) + alpha_xy * C2(i - 1, j + 1, t) + alpha_xy * C2(i + 1, j - 1, t);
        end
    end
    figure(h);
    subplot(121)
    imagesc(C1(:, :, t));
    axis square
    title(['4 neighbor update, ITR = ', num2str(t)]);

    subplot(122)
    imagesc(C2(:, :, t));
    axis square
    title(['8 neighbor update, ITR = ', num2str(t)]);
    
    pause(1);
end

