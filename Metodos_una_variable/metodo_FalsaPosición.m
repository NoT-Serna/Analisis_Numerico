function p = FalsaPosicion(f, p0, p1, TOL, N0)
% FalsaPosicion Encuentra una solución a f(x) = 0 usando el Método de la Falsa Posición.
%   p = FalsaPosicion(f, p0, p1, TOL, N0) encuentra una raíz de la función continua f
%   en el intervalo [p0, p1], donde f(p0) y f(p1) tienen signos opuestos, con una
%   tolerancia TOL y un número máximo de iteraciones N0.

% Paso 1: Inicialización de variables
i = 2;
q0 = f(p0);
q1 = f(p1);

% Paso 2: Bucle principal
while i <= N0
    % Paso 3: Cálculo de la nueva aproximación
    p = p1 - q1 * (p1 - p0) / (q1 - q0);
    
    % Calcula el error relativo si p no es cero
    if p ~= 0
        Er = abs((p - p1) / p);
    else
        Er = abs(p - p1); % Si p es 0, el error relativo es igual a la diferencia absoluta
    end
    
    % Paso 4: Verificación de la condición de parada basada en la diferencia absoluta y el error relativo
    if abs(p - p1) < TOL || Er < TOL
        fprintf('El procedimiento fue exitoso después de %d iteraciones.\n', i-1);
        return
    end
    
    % Paso 5: Actualización del contador de iteraciones
    i = i + 1;
    
    % Paso 5.1: Evaluación de f en la nueva aproximación
    q = f(p);
    
    % Paso 6: Decidir cómo actualizar los puntos p0 y p1
    if q * q1 < 0
        p0 = p1;
        q0 = q1;
    end
