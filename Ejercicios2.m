%----------Pregunta 1------------------%
%Encontrar los intercambios de filas para resolver los sistemas lineales
%usando Gauss con sustitución inversa
% 1. ----------------------
% Definir la matriz de coeficientes y el vector de términos independientes
A = [13, 17, 1; 
    0, 1, 19; 
    0, 12, -1];

b = [5; 1; 0];

x = gauss_eliminacion(A,b);
disp('Solución del sistema:');
disp(x);

clc;
clear;

% 2. ---------------------
A = [1, 1, -1;
    0, 12, -1;
    2, 1, 1;];

b = [0; 4; 5;];

x = gauss_eliminacion(A, b);
disp('Solución del sistema:');
disp(x);

clc;
clear;

%-----Pregunta 2 -------------%
% 1.
A = [2,0,0,0;
    1, 1.5,0,0;
    0,-3,0.5,0;
    2,-2,1,1;];

b = [3; 4.5; -6.6; 0.8];

x = gauss_eliminacion(A, b);
disp('Solución del sistema:');
disp(x);


%------Pregunta 3------------%
%pCO2 = presión parcial del CO2 en la atmósfera calcular el pH de la lluvia
%dado los siguientes parámetros

%Como es un problema no lineal, es necesario usar la matriz jacobiana para
%linearlizar el sistema tal como lo hice en Modeling of Dynamical Systems
clc; clear;

% ---------- CONSTANTES ----------
KH   = 10^-1.46;
K1   = 10^-6.3;
K2   = 10^-10.3;
Kw   = 10^-14;
pCO2 = 375e-6;    % 375 ppm en atm

% ---------- INICIALIZACIÓN ----------
x = [1e-5; 1e-5; 1e-6; 1e-9];  % [H+, HCO3-, CO3--, OH-]
tol_newton = 1e-12;
maxit_newton = 50;

%{
for k = 1:maxit_newton
    H = x(1); HCO3 = x(2); CO3 = x(3); OH = x(4);

    % ---------- SISTEMA NO LINEAL ----------
    F = [
        K1 - 1e6*(H*HCO3)/(KH*pCO2);
        K2 - (H*CO3)/HCO3;
        Kw - H*OH;
        HCO3 + 2*CO3 + OH - H
    ];

    % ---------- JACOBIANA ----------
    J = [
        -1e6 * HCO3 / (KH*pCO2),  -1e6 * H / (KH*pCO2),    0,       0;
        -CO3 / HCO3,               (H*CO3)/(HCO3^2),      -H/HCO3, 0;
        -OH,                       0,                     0,      -H;
        -1,                        1,                     2,       1
    ];

% Pendiente
%}


clc;
clear;

%----Pregunta 4 -------- %
%Ley de la conservación de la masa en tres reactores donde se miden la
%conecntración de algún agente químico reactor 1 = x reactor 2 = y reactor
%3 = z donde se cumple que: ENTRADAS  = SALIDAS

% EC1: 130x - 30y + 0z = 200
% EC2: 90x - 90y + 0z = 0
% EC3: -40 x -60y +120z = 500

A = [130, -30, 0;
    90, -90, 0;
    -40, -60, 120
];

b = [200; 0; 500];

x = gauss_eliminacion(A, b);
disp('Solución del sistema de reactores:');
disp(x);


clc;
clear;

%------Pregunta 5 -----------%
%Calcular el error exacto y el relativo
 A = [58.9, 0.03;
     -6.10, 5.31];

 b = [59.2; 47;];

 x_true = [1; 10];

 x_approx = [1.02; 9.98];

% Error real
error_x = norm(x_true-x_approx, inf);

% Residuo del sistema aproximado
r = b - A*x_approx;

% Norma infinito de A y del residuo
normA_inf = norm(A,inf);
normr_inf = norm(r, inf);

% Condición de A en norma infinito
K_inf = cond(A, inf);

% Error relativo
error_rel = K_inf * (normr_inf / normA_inf);

%Resultados

fprintf('||x-x_approx|| inf = %.6f\n', error_x);
fprintf('K_inf(A) * ||b - A x_tilde||_inf / ||A||_inf = %.6f\n', error_rel);
clc;
clear;

%----------Pregunta 6 ---------------------%
%Dado el sistema lineal

A= [1, 1/2, 1/3;
    1/2, 1/3, 1/4;
    1/3, 1/4, 1/5;];

b = [5/6; 5/12; 17/60;];

x_true = [1; -1; 1;];

%Eliminación Gaussiana con redondeo aritmético a tres cifras:
% Aplicar eliminación gaussiana
x = gauss_eliminacion(A, b);
disp('Solución del sistema con redondeo aritmético:');
disp(x);
%Método del gradiente conjugado
C = diag(diag(A));
n = 100;
tol = 1e-6;
[x_sol, res] = gradiente_conjugado(A, b, C, n, tol);
disp(x_sol)
plot(res, '-o')
xlabel('Iteración')
ylabel('Norma del residuo')
title('Convergencia del método del gradiente conjugado preacondicionado')
grid on

clc;
clear;

%-----Pregunta 7 ----------------%
%Se han reportado ciertos datos que se modelaron a paritr de:
% x = e^(y-b)/a, a y b son parámetros
% 1. Transformación para linearlizar la función
% Transformación logarítmica para linealizar la función
%ln(x) = (y - b) / a;
% y = a ln(x) +b 
% Definir los parámetros a y b a partir de regresión lineal
% Datos
x = [1, 2, 3, 5, 5];
y = [0.5, 2, 2.9, 3.5, 4];

% --- MÉTODO 1: Regresión logarítmica ---
fprintf('=== MÉTODO 1: REGRESIÓN LOGARÍTMICA ===\n');
% Transformación logarítmica
x_lin = log(x);

% Ajuste por regresión lineal
[a0, a1, syx, r2] = regress(x_lin, y);

% Mostrar resultados
fprintf('Recta linealizada: y = %.3f + %.3f ln(x)\n', a0, a1);
fprintf('Error estándar: %.3f\n', syx);
fprintf('Coeficiente de determinación R²: %.4f\n', r2);

% Interpretación con respecto a los parámetros del modelo original
a = a1;    % a = pendiente
fprintf('El parámetro a = %.3f\n', a);
b = a0;    % b = intercepto
fprintf('El parámetro b = %.3f\n', b);

fprintf('\nModelo original: x = exp((y - %.3f)/%.3f)\n', b, a);

%Evaluar y para x = 2.6
x_eval = 2.6;
y_eval_log = a*log(x_eval)+b;
fprintf('Para x = %.2f, y = %.4f\n', x_eval, y_eval_log);

% --- MÉTODO 2: Serie de Fourier truncada USANDO LA FUNCIÓN ---
fprintf('\n=== MÉTODO 2: SERIE DE FOURIER TRUNCADA ===\n');

% Probar diferentes números de términos
for N_terms = 1:2
    fprintf('\n--- Fourier con N = %d términos ---\n', N_terms);
    
    [coeff, y_pred_fourier, r2_fourier, syx_fourier, equation_str] = SerieTruncadaFourier(x, y, N_terms, x_eval);
    
    fprintf('Coeficientes: [');
    fprintf('%.4f ', coeff);
    fprintf(']\n');
    
    fprintf('Ecuación: %s\n', equation_str);
    fprintf('R² = %.4f\n', r2_fourier);
    fprintf('Error estándar = %.4f\n', syx_fourier);
    fprintf('y(%.1f) = %.4f\n', x_eval, y_pred_fourier(end));
    
    % Guardar resultados para N=1 para la comparación visual
    if N_terms == 1
        y_eval_fourier = y_pred_fourier(end);
        fourier_func = @(x_val) coeff(1) + ...
            sum(coeff(2:2:end) .* cos(2*pi*(1:N_terms)'.*(x_val - min(x))/(max(x)-min(x)))) + ...
            sum(coeff(3:2:end) .* sin(2*pi*(1:N_terms)'.*(x_val - min(x))/(max(x)-min(x))));
    end
end

% --- COMPARACIÓN VISUAL ---
figure;

% Subplot 1: Regresión logarítmica
subplot(1,2,1);
x_plot = linspace(min(x), max(x), 100);
y_plot_log = a*log(x_plot) + b;
plot(x, y, 'o', x_plot, y_plot_log, '-r', 'LineWidth', 2);
hold on;
plot(x_eval, y_eval_log, 'sg', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlabel('x'); ylabel('y');
title('Regresión Logarítmica');
legend('Datos', 'Ajuste', 'x=2.6', 'Location', 'best');
grid on;

% Subplot 2: Serie de Fourier
subplot(1,2,2);
y_plot_fourier = arrayfun(fourier_func, x_plot);
plot(x, y, 'o', x_plot, y_plot_fourier, '-b', 'LineWidth', 2);
hold on;
plot(x_eval, y_eval_fourier, 'sg', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlabel('x'); ylabel('y');
title('Serie de Fourier Truncada (N=1)');
legend('Datos', 'Ajuste', 'x=2.6', 'Location', 'best');
grid on;

% --- TABLA COMPARATIVA ---
fprintf('\n=== COMPARACIÓN DE MÉTODOS ===\n');
fprintf('Método            | R²      | Error estándar | y(2.6)\n');
fprintf('----------------------------------------------------\n');
fprintf('Regresión log     | %.4f  | %.4f        | %.4f\n', r2, syx, y_eval_log);
fprintf('Fourier (N=1)     | %.4f  | %.4f        | %.4f\n', r2_fourier, syx_fourier, y_eval_fourier);

clear;
clc;



%---- Pregunta 8 ---------------%
x = [0.5, 1, 2, 3, 4];
y = [10.4, 5.8, 3.3, 2.4, 2];

% Paso 1: Transformación para linearizar
% y = ((a + √x)/(b√x))²
% √y = (a + √x)/(b√x)
% √y = a/(b√x) + 1/b
% Hacemos: z = √y, u = 1/√x
% Obtenemos: z = (a/b)u + 1/b

u = 1./sqrt(x);
z = sqrt(y);

% Paso 2: Aplicar TU función de regresión lineal
[a0, a1, syx, r2] = regress(u, z);

% Paso 3: Encontrar parámetros a y b del modelo original
% De: z = a1*u + a0
% Comparando con: z = (a/b)u + 1/b
% Tenemos: a1 = a/b  y  a0 = 1/b

b = 1 / a0;
a = a1 * b;

% Mostrar resultados
fprintf('=== RESULTADOS DE LA REGRESIÓN ===\n');
fprintf('Parámetros de la regresión lineal:\n');
fprintf('a0 (intercepto) = %.6f\n', a0);
fprintf('a1 (pendiente)  = %.6f\n', a1);
fprintf('Error estándar (syx) = %.6f\n', syx);
fprintf('Coeficiente R² = %.6f\n', r2);

fprintf('\nParámetros del modelo original:\n');
fprintf('a = %.6f\n', a);
fprintf('b = %.6f\n', b);

% Paso 4: Predicción para x = 1.6
x_pred = 1.6;
y_pred = ((a + sqrt(x_pred)) / (b * sqrt(x_pred)))^2;
fprintf('Para x = %.1f, y_pred = %.6f\n', x_pred, y_pred);

clc;
clear;

% ---- Pregunta 9 ------------- %
x1 = [0, 0, 1, 2, 0, 2, 1];
x2 = [0, 2, 2, 4, 4, 6, 2];
y  = [14, 21, 11, 12, 23, 14, 6];

% Matriz de diseño (agregar columna de 1's para el intercepto)
X = [ones(length(x1), 1), x1', x2'];

% Vector de respuestas
Y = y';

% Coeficientes por mínimos cuadrados
B = X \ Y;  % [b0; b1; b2]

% Mostrar resultados
fprintf('Ecuación ajustada: y = %.3f + %.3f*x1 + %.3f*x2\n', B(1), B(2), B(3));

% Predicciones del modelo
y_pred = X * B;

% Calcular R²
SStot = sum((Y - mean(Y)).^2);
SSres = sum((Y - y_pred).^2);
R2 = 1 - SSres / SStot;

fprintf('Coeficiente de determinación R² = %.4f\n', R2);

% Graficar comparación entre datos reales y predichos
figure;
plot(Y, y_pred, 'o', 'MarkerFaceColor', 'b');
xlabel('y real');
ylabel('y predicha');
title('Regresión Lineal Múltiple');
grid on;
axis equal;

clc;
clear;

%-- Problema 10----%
x = [0.2, 0.5, 0.8, 1.2, 1.7, 2, 2.3];
y = [500, 700, 1000, 1200, 2200, 2650, 3750];

%Regresion cuadrática
% Ajuste de regresión cuadrática
X_quad = [ones(length(x), 1), x', x'.^2]; % Matriz de diseño para regresión cuadrática
Y_quad = y'; % Vector de respuestas

% Coeficientes por mínimos cuadrados
B_quad = X_quad \ Y_quad;  % [b0; b1; b2]

% Mostrar resultados
fprintf('Ecuación ajustada: y = %.3f + %.3f*x + %.3f*x^2\n', B_quad(1), B_quad(2), B_quad(3));

% Predicciones del modelo
y_pred_quad = X_quad * B_quad;

% Calcular R²
SStot_quad = sum((Y_quad - mean(Y_quad)).^2);
SSres_quad = sum((Y_quad - y_pred_quad).^2);
R2_quad = 1 - SSres_quad / SStot_quad;

fprintf('Coeficiente de determinación R² = %.4f\n', R2_quad);

%{
% Graficar comparación entre datos reales y predichos
figure;
plot(Y_quad, y_pred_quad, 'o', 'MarkerFaceColor', 'r');
xlabel('y real');
ylabel('y predicha');
title('Regresión Cuadrática');
grid on;
axis equal;
%}

clc;
clear;



%--- Problema 11----- %
%{
% Radiación solar en una ciudad durante el año
% Asumiendo 30 días por mes

% Meses (1 a 12)
meses = 1:12;
radiacion = [144, 188, 245, 311, 351, 359, 308, 287, 260, 211, 159, 131];

% Convertir meses a días (cada mes = 30 días)
dias = (meses - 1) * 30 + 15;  % Día central de cada mes (día 15, 45, 75, ...)

% Ajustar una sinusoide (N=1 término de Fourier)
N_terms = 1;
%[coeff, y_pred, r2, syx, equation_str] = SerieTruncadaFourier(dias, radiacion, N_terms);

fprintf('=== SINUSOIDE AJUSTADA ===\n');
%fprintf('Ecuación: %s\n', equation_str);
fprintf('R² = %.4f\n\n', r2);

% Pronosticar radiación a mediados de agosto (mes 8)
% Agosto: mes 8 → día = (8-1)*30 + 15 = 225
dia_agosto = 225;
radiacion_agosto = coeff(1) + coeff(2)*cos(2*pi*1*(dia_agosto-15)/360) + coeff(3)*sin(2*pi*1*(dia_agosto-15)/360);

fprintf('PRONÓSTICO PARA MEDIADOS DE AGOSTO:\n');
fprintf('Mes: Agosto (mes 8)\n');
fprintf('Día del año: %d\n', dia_agosto);
fprintf('Radiación pronosticada: %.1f W/m²\n', radiacion_agosto);

% Gráfica simple
figure;
meses_nombres = {'E','F','M','A','M','J','J','A','S','O','N','D'};
plot(meses, radiacion, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;

Graficar la sinusoide ajustada

meses_plot = linspace(1, 12, 100);
dias_plot = (meses_plot - 1) * 30 + 15;
radiacion_plot = coeff(1) + coeff(2)*cos(2*pi*1*(dias_plot-15)/360) + coeff(3)*sin(2*pi*1*(dias_plot-15)/360);
plot(meses_plot, radiacion_plot, 'b-', 'LineWidth', 2);

% Marcar el pronóstico de agosto
plot(8, radiacion_agosto, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

xlabel('Mes');
ylabel('Radiación solar (W/m²)');
title('Ajuste Sinusoidal de Radiación Solar');
xticks(1:12);
xticklabels(meses_nombres);
legend('Datos reales', 'Sinusoide ajustada', 'Pronóstico agosto', 'Location', 'best');
grid on;
%}

clc;
clear;

%----Pregunta 12 ----- %
%pH de un reactor varía en formas sinusoidales
% El periodo es de 24h
temp = [0,2,4,5,7,9,12,15,17,19,20,24];
pH   = [7.3, 7, 7.1, 6.5, 7.4, 7.2, 8.2, 8.9, 8.8, 8.9, 8.7, 7.4];
T = 24;

[coef_pH, y_pH, r2_pH, syx_pH, eq_pH] = SerieTruncadaFourier(temp, pH, 1, T);

fprintf("%s\n", eq_pH);
fprintf("R² = %.4f\n", r2_pH);
fprintf("Error estándar = %.4f\n", syx_pH);

clc;
clear;

%----Pregunta 13 ------ %
clear; clc;

a = 1;
b = 2;
h = 0.5;
n = (b - a)/h;

objetivo = log(2);

% Ecuación del sistema
f = @(x, Y) [Y(2); -(Y(2))^2 - Y(1) + log(x)];

% Disparo inicial para y'(1)
s = 0;      
tol = 1e-6;

for iter = 1:20
    
    % Resolver con el disparo actual
    y0 = [0; s];   % y(1)=0
    [t, Y] = runge_kutta_4(f, a, y0, h, n);
    
    % Error
    F = Y(end,1) - objetivo;
    
    if abs(F) < tol
        break;
    end
    
    % Segundo disparo para método secante
    s2 = s + 0.01;
    y0_2 = [0; s2];
    [~, Y2] = runge_kutta_4(f, a, y0_2, h, n);
    
    F2 = Y2(end,1) - objetivo;
    
    % Actualizar con método secante
    s = s - F*(s2 - s)/(F2 - F);
end

disp('Solución aproximada:')
disp(Y(:,1))
