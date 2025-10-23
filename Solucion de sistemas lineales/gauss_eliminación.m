function x = gauss_eliminacion(A, b)
    % A es la matriz de coeficientes
    % b es el vector de términos independientes
    % x es el vector solución del sistema Ax = b

    n = length(b);
    Ab = [A b];  % Matriz aumentada

    % Eliminación hacia adelante
    for k = 1:n-1
        for i = k+1:n
            factor = Ab(i,k) / Ab(k,k);
            Ab(i,k:end) = Ab(i,k:end) - factor * Ab(k,k:end);
        end
    end

    % Sustitución hacia atrás
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (Ab(i,end) - Ab(i,i+1:n) * x(i+1:n)) / Ab(i,i);
    end
end
