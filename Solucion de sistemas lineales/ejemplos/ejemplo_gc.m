
% Definir el sistema A*x = b
A = [4 -1 0; 
    -1 4 -1; 
     0 -1 4];

x = [1; 2; 3];         % Vector de solución verdadera (columna)
b = A * x;             % Calcular el lado derecho correctamente

C = diag(diag(A));     % Preacondicionador simple (diagonal de A)
n = 100;               % Máximo de iteraciones
tol = 1e-6;            % Tolerancia

% Llamar a la función de gradiente conjugado preacondicionado
[x_sol, res] = gradiente_conjugado(A, b, C, n, tol);

% Mostrar la solución obtenida
disp('Solución aproximada:')
disp(x_sol)

% Comparar con la solución verdadera
disp('Solución exacta:')
disp(x)

% Graficar la convergencia
plot(1:length(res), res, '-o', 'LineWidth', 1.5)
xlabel('Iteración')
ylabel('Norma del residuo')
title('Convergencia del método del gradiente conjugado preacondicionado')
grid on
