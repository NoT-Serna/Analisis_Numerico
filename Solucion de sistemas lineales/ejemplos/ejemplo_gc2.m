%Desarollar la función de energía potencial  para el sistema, Minimizar la
%función de energía para determinar los desplazamientos en el equilibrio x1
%y x2 F = 100N Ka=20N/m Kb = 15N/m

% Ecuaciones Lineales
% 1. (Ka+Kb)x1-Kbx2=0
% 2. -Kbx1+Kbx2 = F
% Parámetros
Ka = 20;
Kb = 15;
F = 100;


A = [Ka + Kb, -Kb; 
    -Kb, Kb];
b = [0; F];
x = A \ b;

x1 = x(1);
x2 = x(2);

fprintf('x1 = %.4f m\n', x1)
fprintf('x2 = %.4f m\n', x2)

% ---- Energía potencial total en equilibrio ----
U = 0.5*Ka*x1^2 + 0.5*Kb*(x2 - x1)^2 - F*x2;
fprintf('Energía potencial mínima = %.4f J\n', U)
