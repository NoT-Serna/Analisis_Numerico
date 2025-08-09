function p = NewtonRaphson(f, df, p0, TOL, N0)
    %Inicializacion
    i = 1;
    p = p0;

    % Bucle principal
    while i <= N0
        % Calcula la nueva aproximacion 
        p_new = p -f(p)/df(p);

        % Calcula el error relativo
        if p_new ~= 0
            Er = abs((p_new -p)/p_new);
        else
            % Evita la division por cero si p_new es cero
            Er = abs(p_new -p);
        end

        % Verifica si el error relativo es menor que a tolerancia 

        if Er < TOL
            fprintf(['La raíz aproximada es: %.5f después de %d iteraciones.' ...
                '\n'], p_new, i)
            return
        end
        %Actualiza i para la siguiente iteracion
        i = i + 1;
        p = p_new;
    end
    % Si se alcanza el número máximo de iteraciones sin converger
    fprintf('El método falló después de %d iteraciones.\n', N0);
end
