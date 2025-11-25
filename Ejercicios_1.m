%----- Problema 1 -------%
%Modelo Michaelis-Menten para enzimas alostéricas donde:
%Datos:
S = [1.3 1.8 3 4.5 6 8 9];
v = [0.07 0.13 0.22 0.275 0.335 0.35 0.36];

% Linealización
X = 1./(S.^2);
Y = 1./v;

% Matriz del sistema
A = [X' ones(length(S),1)];

% SVD
[U,Sigma,V] = svd(A,'econ');
c = V * diag(1./diag(Sigma)) * U' * Y';

a = c(1);
b = c(2);

% Parámetros del modelo
Vmax = 1/b;
Ks = sqrt(a * Vmax);

fprintf("Vmax = %.6f\n", Vmax);
fprintf("Ks   = %.6f\n", Ks);

% Predicción del modelo
v_pred = (Vmax .* S.^2) ./ (Ks^2 + S.^2);

% ---- Cálculo de R² ----
SS_res = sum((v - v_pred).^2);
SS_tot = sum((v - mean(v)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf("R² = %.6f\n", R2);

% ---- Gráfica ----
figure;
plot(S, v, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
S_plot = linspace(min(S), max(S), 200);
v_curve = (Vmax .* S_plot.^2) ./ (Ks^2 + S_plot.^2);
plot(S_plot, v_curve, 'r-', 'LineWidth', 2);
grid on;

xlabel('S');
ylabel('v');
title(sprintf('Ajuste Michaelis–Menten (R² = %.4f)', R2));
legend('Datos experimentales','Curva ajustada','Location','SouthEast');

clear;
clc;

%----- Problema 2 -------%
% Promedios del Dow Jones (DIJA) entre marzo del 2013 y junio 2014
x = 0:1:15;
y = [14090, 14573, 14701, 15254, 14975, 15628, 14834, 15193,...
     15616, 16009, 16441, 15373, 16168, 16533, 16559, 16744];

% (a) Ajuste polinómico grado 4
p = polyfit(x, y, 4);
y_fit = polyval(p, x);

% (b) Aproximaciones para abril 2013 (x = 8) y abril 2014 (x = 13)
y_8  = polyval(p, 8);
y_13 = polyval(p, 13);

fprintf("Aproximación para abril 2013 (x=8): %.2f\n", y_8);
fprintf("Aproximación para abril 2014 (x=13): %.2f\n\n", y_13);

% (c) Comparación con valores reales
real_2013 = 14613;
real_2014 = 16256;

error_2013 = abs((y_8 - real_2013) / real_2013) * 100;
error_2014 = abs((y_13 - real_2014) / real_2014) * 100;

fprintf("Valor real 2013: %d - Aproximado: %.2f - Error: %.2f%%\n", ...
        real_2013, y_8, error_2013);

fprintf("Valor real 2014: %d - Aproximado: %.2f - Error: %.2f%%\n\n", ...
        real_2014, y_13, error_2014);

% (d) Predicción para 17 junio 2014 (aprox. x = 15.5)
y_15_5 = polyval(p, 15.5);
fprintf("Aproximación para el 17 de junio 2014 (x=15.5): %.2f\n", y_15_5);

clc;
clear;


%----Problema 3 -------%
% Problema 3 - Shooting Lineal
f = @(r, u) [u(2); -(2/r)*u(2)];

R1 = 2;
R2 = 4;
V1 = 110;
V2 = 0;

% Vector de r fijo para ambas integraciones
r_span = linspace(R1, R2, 100);

% Disparo A: uA(R1) = V1, uA'(R1) = 0
u0_a = [V1; 0];
[~, ua] = ode45(f, r_span, u0_a);

% Disparo B: uB(R1) = 0, uB'(R1) = 1
u0_b = [0; 1];
[~, ub] = ode45(f, r_span, u0_b);

% Calcular alpha
alpha = (V2 - ua(end,1)) / ub(end,1);

% Evaluar en r = 3
u3_approx = interp1(r_span, ua(:,1), 3) + alpha * interp1(r_span, ub(:,1), 3);

disp(['Aproximación u(3) = ', num2str(u3_approx)]);

% Solución exacta para comparar
u_exact = @(r) (V1 * R1 ./ r) .* ((R2 - r) / (R2 - R1));
disp(['Valor exacto u(3) = ', num2str(u_exact(3))]);
disp(['Error = ', num2str(abs(u3_approx - u_exact(3)))]);
