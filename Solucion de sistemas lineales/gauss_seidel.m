function [x, iter] = gauss_seidel(A, b, x0, tol, max_iter)
    % Entrada:
    % A - Matriz de coeficientes (n x n)
    % b - Vector de términos independientes (n x 1)
    % x0 - Aproximación inicial (n x 1)
    % tol - Tolerancia para la convergencia
    % max_iter - Máximo número de iteraciones permitidas

    % Salida:
    % x - Solución aproximada (n x 1)
    % iter - Número de iteraciones realizadas

    n = length(b); % Número de ecuaciones (o tamaño de la matriz A)
    x = x0; % Inicializamos x con la aproximación inicial
    iter = 0; % Inicializamos el contador de iteraciones

    for k = 1:max_iter
        x_old = x; % Guardamos la solución anterior para comparar después
        
        % Iterar sobre cada ecuación
        for i = 1:n
            sum1 = 0;
            sum2 = 0;
            
            % Primera sumatoria: j de 1 hasta i-1
            for j = 1:i-1
                sum1 = sum1 + A(i,j) * x(j);
            end
            
            % Segunda sumatoria: j de i+1 hasta n
            for j = i+1:n
                sum2 = sum2 + A(i,j) * x0(j); % Usamos x0 para las iteraciones anteriores
            end
            
            % Actualización de x_i con la fórmula de Gauss-Seidel
            x(i) = (b(i) - sum1 - sum2) / A(i,i);
        end
        
        % Comprobar la convergencia
        if norm(x - x_old, inf) < tol
            disp('El procedimiento fue exitoso.');
            return;
        end
        
        % Actualizar la aproximación inicial para la siguiente iteración
        x0 = x;
        iter = iter + 1;
    end
    
    % Si se excede el número máximo de iteraciones
    disp('Se ha excedido el número máximo de iteraciones.');
end
