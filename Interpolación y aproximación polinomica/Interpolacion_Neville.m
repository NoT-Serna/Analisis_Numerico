classdef NevilleInterp
    %NevilleInterp Realiza interpolación polinómica usando el algoritmo de Neville.
    %
    %   Esta clase implementa el algoritmo de Neville como se describe en el
    %   libro "Numerical Recipes". Permite encontrar el valor de un polinomio
    %   que pasa a través de un conjunto de puntos de datos.
    %
    %   La clase puede operar de dos maneras:
    %   1.  **Interpolación Global:** Si se usa el número total de puntos,
    %       calcula el único polinomio de grado N-1 que pasa por los N puntos.
    %   2.  **Interpolación Local:** Si se especifica un número M < N de puntos,
    %       realiza una interpolación local usando los M puntos más cercanos
    %       al punto de evaluación.
    %
    %   Uso Básico:
    %       % 1. Crear un objeto con los datos (x, y) y el orden (M)
    %       interpolador = NevilleInterp(x_datos, y_datos, M);
    %
    %       % 2. Interpolar un valor en un punto x
    %       [y_valor, interpolador] = interpolador.interp(x);
    %
    %       % 3. Obtener la estimación del error del último cálculo
    %       error_estimado = interpolador.dy;

    % --- PROPIEDADES ---
    % Las propiedades almacenan el estado del objeto interpolador.
    properties
        xx (:,1) double % Vector columna con las coordenadas 'x' de los datos.
        yy (:,1) double % Vector columna con las coordenadas 'y' de los datos.
        m  (1,1) double % Número de puntos a usar en la interpolación (orden + 1).
        dy (1,1) double % Almacena la estimación del error del último valor interpolado.
    end

    % --- MÉTODOS ---
    methods
        function obj = NevilleInterp(xv, yv, m_in)
            %NEVILLEINTERP Constructor de la clase.
            %
            %   Sintaxis:
            %       obj = NevilleInterp(xv, yv, m_in)
            %
            %   Entradas:
            %       xv   - Vector (fila o columna) con las coordenadas x.
            %       yv   - Vector (fila o columna) con las coordenadas y.
            %       m_in - Entero que especifica el número de puntos a usar
            %              en cada interpolación (grado del polinomio + 1).
            %
            %   Salida:
            %       obj - Una nueva instancia de la clase NevilleInterp.

            % Validación de las entradas para evitar errores comunes.
            if numel(xv) ~= numel(yv)
                error('Los vectores de entrada x e y deben tener el mismo tamaño.');
            end
            if m_in > numel(xv) || m_in < 2
                error('M debe ser al menos 2 y no mayor que el número total de puntos.');
            end

            % Asignación de valores a las propiedades del objeto.
            % Se asegura que los vectores sean columnas con `(:)`.
            obj.xx = xv(:);
            obj.yy = yv(:);
            obj.m = m_in;
            obj.dy = 0.0; % Inicialización del error estimado.
        end

        function [y, obj] = interp(obj, x)
            %INTERP Interpola un valor 'y' en un punto 'x'.
            %
            %   Este es el método principal que el usuario llama. Determina
            %   qué subconjunto de puntos usar y luego invoca a `rawinterp`
            %   para realizar el cálculo.
            %
            %   Entradas:
            %       obj - El objeto de la clase.
            %       x   - El punto (escalar) en el que se desea interpolar.
            %
            %   Salidas:
            %       y   - El valor interpolado en el punto x.
            %       obj - El objeto actualizado (importante porque 'dy' cambia).

            % Si 'm' es igual al total de puntos, se usa un polinomio global.
            if obj.m == numel(obj.xx)
                jl = 1; % El subconjunto de datos comienza en el primer índice.
            else
                % Si no, se hace una interpolación local.
                % 1. Encontrar el índice del punto más cercano a 'x'.
                [~, j_closest] = min(abs(obj.xx - x));
                
                % 2. Calcular el índice inicial 'jl' del subconjunto de 'm' puntos,
                %    intentando mantener 'x' centrado.
                jl = j_closest - floor((obj.m - 1) / 2);

                % 3. Asegurarse de que el subconjunto no se salga de los límites.
                if jl < 1
                    jl = 1;
                end
                if jl + obj.m - 1 > numel(obj.xx)
                    jl = numel(obj.xx) - obj.m + 1;
                end
            end
            % Llamar al método interno que implementa el algoritmo.
            [y, obj] = obj.rawinterp(jl, x);
        end

        function [y, obj] = rawinterp(obj, jl, x)
            %RAWINTERP Implementación central del algoritmo de Neville.
            %   (Método interno, no diseñado para ser llamado por el usuario).
            %
            %   Este método construye el "cuadro" de Neville para un
            %   subconjunto específico de datos y calcula el valor interpolado.

            % Extraer el subconjunto de datos (x e y) a utilizar.
            xa = obj.xx(jl : jl + obj.m - 1);
            ya = obj.yy(jl : jl + obj.m - 1);

            % 'c' y 'd' son los vectores que almacenan las correcciones en cada
            % Se inicializan con los valores de 'y'.
            c = ya;
            d = ya;
            
            % Encontrar el punto más cercano dentro del subconjunto.
            [~, ns_matlab] = min(abs(x - xa));
            
            % La primera aproximación al resultado es el valor 'y' del punto más cercano.
            y = ya(ns_matlab);
            
            % Se prepara 'ns_c' (índice 0-based, estilo C) para el bucle.
            % La línea `y=ya[ns--]` en C++ usa 'ns' y LUEGO decremento.
            % Aquí se emula ese comportamiento para la primera iteración del bucle.
            ns_c = ns_matlab - 1; % Convertir a 0-based.
            ns_c = ns_c - 1;      % Emular el post-decremento.
            
            % Bucle principal: construye el cuadro de Neville columna a columna.
            for m_loop = 1:(obj.m - 1)
                % Bucle interno: calcula cada elemento de la columna actual.
                for i = 1:(obj.m - m_loop)
                    % Diferencias entre los puntos x y el punto de evaluación.
                    ho = xa(i) - x;
                    hp = xa(i + m_loop) - x;

                    % 'w' es la diferencia entre los valores de la columna anterior.
                    w = c(i + 1) - d(i);
                    den = ho - hp;
                    
                    % Evitar división por cero si dos puntos 'x' son idénticos.
                    if den == 0.0
                        error('Error en NevilleInterp: puntos xa idénticos.');
                    end

                    % Calcular y actualizar los valores de 'c' y 'd' para esta
                    % columna, que son las correcciones de orden superior[cite: 49, 50, 51].
                    den = w / den;
                    d(i) = hp * den;
                    c(i) = ho * den;
                end

                % Decidir qué corrección añadir a 'y'. Se elige la ruta a través
                % del cuadro que mantiene el resultado centrado en 'x'.
                % Esta es la traducción de la línea:
                % y += (dy = (2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--])); 
                if ( 2*(ns_c + 1) < (obj.m - m_loop) )
                    % Ruta "hacia abajo" en el cuadro: usa 'c'.
                    % El índice en MATLAB es (ns_c+1)+1 = ns_c+2.
                    obj.dy = c(ns_c + 2);
                else
                    % Ruta "hacia arriba" en el cuadro: usa 'd'.
                    % El índice en MATLAB es ns_c+1.
                    obj.dy = d(ns_c + 1);
                    % Se decrementa el índice para la siguiente iteración,
                    % emulando el `d[ns--]` del C++.
                    ns_c = ns_c - 1;
                end
                % Se añade la corrección de este nivel al resultado acumulado.
                y = y + obj.dy;
            end
        end
    end
end

% ------------------------------------------USO DEL ALGORITMO-------------------------------------
% x = [2,4,1,6,7,9]; y = [2,1,5,7,3,8]; 
% interp = NevilleInterp(x,y,4); El valor m_in representa el grado del polinomio y los puntos a usar en la interpolacion si m = # numero de puntos interpolación global si m < # num puntos es interpolacion local
% [val, interp] = interp.interp(12.4); Ese 12.4 es el valor a interpolar
% disp(val); valor al que correponde la interpolacion del punto definido como 12.4, la interpolacion va tener mucho error si el punto esta afuera del rango de puntos, si se usan muchos puntos se puede causar el fenomenode Runge donde tambien se peride percisión por "sobreestimar o sobreinterpolar puntos"
% disp(interp.dy) ultima correcion añadida/estimación , usada para caluclar el error realtivo dy/val en abs val, ese valor represneta el error de la estimacion con respecto al val, mientras mas grande peor la estimacion del valor interpolador en ese punto
