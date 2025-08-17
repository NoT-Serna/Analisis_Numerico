function p = FalsaPosicion(f, p0, p1, TOL, N0)
% FalsaPosicion Encuentra una solución a f(x) = 0 usando el Método de la Falsa Posición.
%   p = FalsaPosicion(f, p0, p1, TOL, N0) encuentra una raíz de la función continua f
%   en el intervalo [p0, p1], donde f(p0) y f(p1) tienen signos opuestos, con una
%   tolerancia TOL y un número máximo de iteraciones N0.

    % Paso 1: Inicialización
    i = 2;
    q0 = f(p0);
    q1 = f(p1);

    % Paso 2: Bucle principal
    while i <= N0
        % Paso 3: Cálculo de la nueva aproximación (interpolación lineal)
        p = p1 - q1 * (p1 - p0) / (q1 - q0);

        % Calcula el error relativo si p no es cero
        if p ~= 0
            Er = abs((p - p1) / p);
        else
            Er = abs(p - p1);
        end

        % Paso 4: Condición de parada
        if abs(p - p1) < TOL || Er < TOL
            fprintf('La raíz aproximada es: %.5f\n', p)
            fprintf('El procedimiento fue exitoso después de %d iteraciones.\n', i-1)
            return
        end

        % Paso 5: Evaluar f en la nueva aproximación
        q = f(p);

        % Paso 6: Actualizar el intervalo
        if q * q0 < 0
            p1 = p;
            q1 = q;
        else
            p0 = p;
            q0 = q;
        end

        % Paso 7: Incrementar contador
        i = i + 1;
    end

    % Si no converge en N0 pasos
    error('El método falló después de %d iteraciones.', N0)
end
