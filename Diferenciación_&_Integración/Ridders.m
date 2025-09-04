function [derivative, error_estimate] = dfridr(func, x, h)
%DFRIDR Computes the derivative of a function using Ridders' method.
%
%   Sintaxis:
%       [deriv, err] = dfridr(func, x, h)
%
%   Descripción:
%       Calcula la derivada de la función 'func' en el punto 'x' usando
%       el método de extrapolación polinomial de Ridders.
%
%   Entradas:
%       func - Handle a la función a diferenciar. Ej: @(x) sin(x).
%       x    - El punto en el cual se evaluará la derivada.
%       h    - Un tamaño de paso inicial. No necesita ser muy pequeño.
%
%   Salidas:
%       derivative     - La estimación de la derivada.
%       error_estimate - Una estimación del error absoluto de la derivada.

% --- Validación de la entrada ---
if h == 0
    error('El tamaño de paso inicial h no puede ser cero.');
end

% --- Constantes del algoritmo ---
NTAB = 10;          % Número máximo de iteraciones (filas en la tabla).
CON = 1.4;          % Factor para reducir el paso h en cada iteración.
CON2 = CON^2;       % Cuadrado del factor, usado en la fórmula de extrapolación.
SAFE = 2.0;         % Factor de seguridad para la condición de salida.

% --- Inicialización ---
a = zeros(NTAB, NTAB); % Pre-alocamos la tabla de extrapolación (tableau).
hh = h;
error_estimate = inf; % Inicializamos el error con un valor grande.

% --- Primera estimación (la menos precisa) ---
% Se calcula la diferencia central con el paso inicial hh.
% Corresponde a a[0][0] en el código C++.
a(1, 1) = (func(x + hh) - func(x - hh)) / (2 * hh);

% --- Bucle principal para construir la tabla de extrapolación ---
for i = 2:NTAB
    % Reducir el tamaño del paso.
    hh = hh / CON;
    
    % Calcular una nueva aproximación con el paso más pequeño.
    % Esto llena la primera fila de la tabla: a(1, i).
    % Corresponde a a[0][i] en C++.
    a(1, i) = (func(x + hh) - func(x - hh)) / (2 * hh);
    
    fac = CON2;
    
    % Bucle de extrapolación (columnas de la tabla).
    % Usa las aproximaciones anteriores para calcular unas nuevas y más precisas.
    for j = 2:i
        % Fórmula de extrapolación de Neville.
        % Corresponde a a[j][i] = (a[j-1][i]*fac - a[j-1][i-1]) / (fac-1.0);
        a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1);
        fac = fac * CON2;
        
        % Estimar el error de esta nueva extrapolación.
        % Se compara la nueva aproximación con las dos que la generaron.
        errt = max(abs(a(j, i) - a(j-1, i)), abs(a(j, i) - a(j-1, i-1)));
        
        % Si el error ha disminuido, guardamos esta como la mejor respuesta.
        if errt <= error_estimate
            error_estimate = errt;
            derivative = a(j, i);
        end
    end
    
    % Condición de salida: Si el error de las aproximaciones de mayor orden
    % (la diagonal de la tabla) empieza a aumentar, significa que los errores
    % de redondeo están dominando. Salimos para evitar un resultado peor.
    if abs(a(i, i) - a(i-1, i-1)) >= SAFE * error_estimate
        return; % Sale de la función y devuelve el mejor resultado encontrado.
    end
end

end

%{
----- USO---------------

f = @(x) -0.1*x.^4 - 0.15*x.^3 - 0.5*x.^2-0.25*x+1.2

f =

  function_handle with value:

    @(x)-0.1*x.^4-0.15*x.^3-0.5*x.^2-0.25*x+1.2

>> [derivative, error_estimate] = dfridr(f,0.5,0.5);
>> disp(derivative)
   -0.9125

>> disp(error_estimate)
%{
