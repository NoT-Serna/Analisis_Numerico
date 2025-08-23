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

%{
-------------- USO----------------------------------------
% Your data points
xData = [-1 0 1];
yData = [1 0 1];

% Create interpolator object
interp = LagrangeInterp(xData, yData);

% Evaluate the interpolating polynomial at some fine grid
xx = linspace(-1,1,200);
yy = interp.lagrange(xx);

% Plot
figure;
plot(xData, yData, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % original points
hold on;
plot(xx, yy, 'b-', 'LineWidth', 1.5); % interpolating polynomial
grid on;
title('Lagrange Interpolation');
legend('Data points','Interpolating polynomial');

%}
