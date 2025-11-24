function [coeff, y_pred, r2, syx, equation_str] = SerieTruncadaFourier(x, y, N_terms, x_eval)
% SERIETRUNCADAFOURIER - Ajusta una serie de Fourier truncada a datos

    % Validación de inputs
    if length(x) ~= length(y)
        error('Los vectores x e y deben tener la misma longitud');
    end
    
    if N_terms < 1 || floor(N_terms) ~= N_terms
        error('N_terms debe ser un entero positivo');
    end
    
    % Definir el intervalo para la serie de Fourier
    L1 = min(x);
    L2 = max(x);
    T = L2 - L1;  % Periodo
    
    % Construir matriz del sistema
    A = zeros(length(x), 2*N_terms + 1);
    
    % Término constante
    A(:,1) = ones(length(x),1);
    
    % Términos coseno y seno
    for n = 1:N_terms
        A(:,2*n) = cos(2*pi*n*(x - L1)/T);
        A(:,2*n+1) = sin(2*pi*n*(x - L1)/T);
    end
    
    % Verificar si el sistema está sobredeterminado
    if size(A,1) < size(A,2)
        warning('Sistema subdeterminado: más parámetros que puntos de datos');
    end
    
    % Resolver por mínimos cuadrados (usar pseudoinversa para estabilidad)
    coeff = pinv(A) * y';
    
    % Calcular predicciones en puntos originales
    y_pred_original = A * coeff;
    
    % Calcular métricas de error
    SSE = sum((y - y_pred_original').^2);
    SST = sum((y - mean(y)).^2);
    r2 = 1 - SSE/SST;
    syx = sqrt(SSE/(length(x) - (2*N_terms + 1)));
    
    % Construir string de la ecuación
    equation_str = sprintf('y = %.4f', coeff(1));
    for n = 1:N_terms
        cos_coeff = coeff(2*n);
        sin_coeff = coeff(2*n+1);
        
        % Formatear término coseno
        if cos_coeff >= 0
            equation_str = sprintf('%s + %.4f cos(2π*%d*(x-%.2f)/%.2f)', ...
                equation_str, cos_coeff, n, L1, T);
        else
            equation_str = sprintf('%s - %.4f cos(2π*%d*(x-%.2f)/%.2f)', ...
                equation_str, abs(cos_coeff), n, L1, T);
        end
        
        % Formatear término seno
        if sin_coeff >= 0
            equation_str = sprintf('%s + %.4f sin(2π*%d*(x-%.2f)/%.2f)', ...
                equation_str, sin_coeff, n, L1, T);
        else
            equation_str = sprintf('%s - %.4f sin(2π*%d*(x-%.2f)/%.2f)', ...
                equation_str, abs(sin_coeff), n, L1, T);
        end
    end
    
    % Evaluar en puntos adicionales si se solicitan
    if nargin > 3
        % Función de evaluación
        fourier_func = @(x_val) coeff(1) + ...
            sum(coeff(2:2:end) .* cos(2*pi*(1:N_terms)'.*(x_val - L1)/T)) + ...
            sum(coeff(3:2:end) .* sin(2*pi*(1:N_terms)'.*(x_val - L1)/T));
        
        % Evaluar en puntos originales y el punto adicional
        y_eval_points = arrayfun(fourier_func, [x, x_eval]);
        y_pred = y_eval_points;
    else
        y_pred = y_pred_original';
    end
end
