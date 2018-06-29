% Test the implementation of the mofified van Der Pol equation
clear
close all;

tspan = [0 800];
% y0 = [.7424051892 0]';
y0 = [2 0]';
[t,y] = ode45('ModifiedVanDerPol', tspan, y0);
% % % [t,y] = ode23('ModifiedVanDerPol', tspan, y0);

figure;
plot(y(:,1),y(:,2));
grid;

figure;
plot(t,y(:,1))
xlabel('n')
ylabel('solution y')
title('Modified van der Pol Equation');
grid

