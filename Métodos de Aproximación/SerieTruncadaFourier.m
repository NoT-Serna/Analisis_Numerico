function [coef, y_pred, R2, syx, eq_str] = SerieTruncadaFourier(x, y, N, T, x_eval)
% Ajuste de serie truncada de Fourier con periodo T (parámetro)
% x : valores de tiempo
% y : datos observados
% N : número de términos de seno/coseno
% T : periodo
% x_eval (opcional) : punto donde evaluar la función

w = 2*pi/T;

x = x(:);
y = y(:);
n_data = length(x);

% Construcción de la matriz A del ajuste
A = ones(n_data, 1);

for n = 1:N
    A = [A  cos(n*w*x)  sin(n*w*x)];
end

% Resolver por mínimos cuadrados
coef = A \ y;

% Predicción en los puntos originales
y_pred = A * coef;

% Estadísticos
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res / SS_tot;
syx = sqrt(SS_res / (n_data - length(coef)));

% Construir ecuación en formato texto
eq_str = sprintf('y = %.4f', coef(1));
k = 2;

for n = 1:N
    eq_str = sprintf('%s + %.4f*cos(%d*(2π/%d)*t) + %.4f*sin(%d*(2π/%d)*t)', ...
        eq_str, coef(k), n, T, coef(k+1), n, T);
    k = k + 2;
end

% Evaluación opcional
if nargin == 5
    A_eval = 1;
    for n = 1:N
        A_eval = [A_eval cos(n*w*x_eval) sin(n*w*x_eval)];
    end
    y_pred = [y_pred; A_eval * coef];
end

end
