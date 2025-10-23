function [x, iter] = jacobi(A, b, x0, tol, max_iter)
    % Entrada:
    % A - Matriz de coeficientes (n x n)
    % b - Vector de términos independientes (n x 1)
    % x0 - Aproximación inicial (n x 1)
    % tol - Tolerancia para la convergencia
    % max_iter - Máximo número de iteraciones
    
    % Salida:
    % x - Solución aproximada (n x 1)
    % iter - Número de iteraciones realizadas

    n = length(b); % Número de ecuaciones (o tamaño de la matriz A)
    x = x0; % Inicializamos x con la aproximación inicial
    iter = 0; % Inicializamos el contador de iteraciones
    for k = 1:max_iter
        x_new = zeros(n, 1); % Crear un nuevo vector para la iteración actual
        for i = 1:n
            sum_ = 0;
            for j = 1:n
                if j ~= i
                    sum_ = sum_ + A(i, j) * x0(j);
                end
            end
            x_new(i) = (b(i) - sum_) / A(i, i);
        end
        
        % Comprobar la convergencia
        if norm(x_new - x0, inf) < tol
            x = x_new;
            disp('El procedimiento fue exitoso.');
            return;
        end
        
        % Actualizar la aproximación para la siguiente iteración
        x0 = x_new;
        iter = iter + 1;
    end
    
    % Si el número máximo de iteraciones fue excedido
    x = x_new;
    disp('Se ha excedido el número máximo de iteraciones.');
end


%{
  USO 
 disp('Solución aproximada:');
disp(x_approx);

disp('Número de iteraciones:');
disp(iter);

% Compare with the true solution
error = norm(x_approx - x_true, inf);
fprintf('Error infinito: %.6e\n', error);


%}
