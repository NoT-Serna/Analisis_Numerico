function p = biseccion(f, a, b, TOL, N0)
    %inicializacion
    i = 1;
    FA = f(a);
    p = a; % Inicializar para evitar division por cero

    % Bucle del metodo de biseccion
    while i <= N0
        % Paso previo: Guardar el valor anterior de p
        p_old = p;

        % Paso 3: Calculo del punto medio
        p = a + (b -a)/2;
        FP = f(p);
        
        % Verificar si p es cero para evitar la division por cero
        if p ~= 0
            Er = abs((p - p_old)/p);
        else
            % Si p es cero use la diferencia absoluta
            Er = abs(b -a);
        end

        % Paso 4: Verificar la condicion de parada
        if FP == 0 || Er < TOL
            fprintf('La raíz aproximada es: %.5f\n', p)
            return
        end
        % Paso 5: Incrementar el contador de iteraciones
        i = i + 1;

        % Paso 6: Actualizar los valores de a o b
        if FA * FP > 0
            a = p;
            FA = FP;
        else 
            b = p;
        end
    end
    % Paso 7: Si se alcanza el numero maximo de iteraciones
    error('El método falló despues del maximo de iteraciones')
end
