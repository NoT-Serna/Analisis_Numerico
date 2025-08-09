function P = lagrangeCoefficients(x, y)
% lagrangeCoefficients Calcula los coeficientes del polinomio de interpolación de Lagrange.
%
%   SINOPSIS:
%       P = lagrangeCoefficients(x, y)
%
%   DESCRIPCIÓN:
%       Esta función calcula los coeficientes del polinomio P(x) que pasa
%       exactamente a través de un conjunto de puntos de datos (x, y). El
%       polinomio resultante tiene un grado de n-1, donde n es el número
%       de puntos.
%
%       El método se basa en la construcción de polinomios base de Lagrange L_i(x)
%       tal que:
%       P(x) = y_1*L_1(x) + y_2*L_2(x) + ... + y_n*L_n(x)
%
%   PARÁMETROS DE ENTRADA:
%       x - Vector fila (1xn) con las coordenadas x de los puntos.
%       y - Vector fila (1xn) con las coordenadas y de los puntos.
%           Debe tener el mismo número de elementos que x.
%
%   PARÁMETROS DE SALIDA:
%       P - Vector fila (1xn) que contiene los coeficientes del polinomio
%           interpolador, ordenados de la potencia más alta a la más baja.
%           Se puede usar con la función polyval.
%
%   EJEMPLO DE USO:
%       % Puntos de datos
%       x = [1, 2, 4];
%       y = [3, -1, 5];
%
%       % Calcular los coeficientes del polinomio
%       P = lagrangeCoefficients(x, y);
%
%       % P debería ser [2.5000, -11.5000, 12.0000]
%       % Esto representa el polinomio 2.5*x^2 - 11.5*x + 12
%
%       % Para verificar, evaluamos el polinomio en los puntos x
%       y_eval = polyval(P, x); % Debería dar [3, -1, 5]
%
%   Ver también: polyfit, polyval.

    % --- Validación de Entradas ---
    if length(x) ~= length(y)
        error('Los vectores de entrada x e y deben tener la misma longitud.');
    end
    
    n = length(x); % Número de puntos de datos
    
    % --- Prealocación de la Matriz de Coeficientes Base ---
    % Cada fila i de L contendrá los coeficientes del polinomio base L_i(x)
    % El polinomio tiene grado n-1, por lo que necesita n coeficientes.
    L = zeros(n, n);

    % --- Construcción de la Matriz de Polinomios Base de Lagrange ---
    % Itera sobre cada punto para construir su polinomio base L_i(x)
    for i = 1:n
        % Inicializa el polinomio base para L_i(x)
        % Empezamos con el polinomio '1' para poder usar conv()
        Li = 1;
        
        % Construye el numerador del polinomio base L_i(x)
        for j = 1:n
            if i ~= j
                % Multiplica por el término (x - x_j)
                % El polinomio [1, -x(j)] representa (x - x_j)
                Li = conv(Li, [1, -x(j)]);
                
                % Divide por el denominador (x_i - x_j)
                Li = Li / (x(i) - x(j));
            end
        end
        
        % Almacena los coeficientes de L_i(x) en la fila i de la matriz L
        L(i, :) = Li;
    end
    
    % --- Cálculo de los Coeficientes del Polinomio Final ---
    % P(x) = sum_{i=1 to n} y_i * L_i(x)
    % Esto se logra multiplicando el vector y por la matriz de coeficientes L
    P = y * L;
end
