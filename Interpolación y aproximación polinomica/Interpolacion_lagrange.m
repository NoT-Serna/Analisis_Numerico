classdef LagrangeInterp < handle
    %LAGRANGEINTERP Implementa interpolación polinómica global
    %mediante el método de Lagrange.
    %
    %   obj = LagrangeInterp(xData, yData)
    %       Construye un objeto con nodos (xData, yData).
    %
    %   yVals = obj.lagrange(xVals)
    %       Devuelve los valores del polinomio de Lagrange que interpola
    %       TODOS los nodos, evaluados en xVals (escalares o vector).
    
    properties
        xv  % nodos x
        yv  % nodos y
    end
    
    methods
        %----------------------------------------------------------
        % Constructor
        %----------------------------------------------------------
        function obj = LagrangeInterp(xData, yData)
            if length(xData) ~= length(yData)
                error('LagrangeInterp: xData e yData deben tener igual longitud.');
            end
            % guardamos como vectores-fila
            obj.xv = xData(:).';
            obj.yv = yData(:).';
        end
        
        %----------------------------------------------------------
        % Método principal de Lagrange
        %----------------------------------------------------------
        function yVals = lagrange(obj, xVals)
            % LAGRANGE Evalúa el polinomio interpolador de Lagrange
            % en los puntos xVals.
            %
            %   yVals = obj.lagrange(xVals)
            %
            %   - xVals : escalar o vector de puntos donde evaluar
            %   - yVals : vector con los valores interpolados
            
            n = length(obj.xv);
            xData = obj.xv;
            yData = obj.yv;
            
            % Prealocación para la salida
            yVals = zeros(size(xVals));
            
            % Evaluamos en cada x solicitado
            for kx = 1:numel(xVals)
                x = xVals(kx);
                P = 0; % Inicializamos el polinomio
                
                % Construcción del polinomio de Lagrange
                for i = 1:n
                    % Calculamos el coeficiente L_i(x)
                    L = 1;
                    for j = 1:n
                        if j ~= i
                            L = L * (x - xData(j)) / (xData(i) - xData(j));
                        end
                    end
                    % Sumamos el término correspondiente
                    P = P + yData(i) * L;
                end
                
                % Guardamos el resultado
                yVals(kx) = P;
            end
        end
    end
end
