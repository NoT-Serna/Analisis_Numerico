function splineFunc = createCubicSpline(x, y)
    % Calcula los coeficientes del spline cúbico
    [a, b, c, d] = cubicSpline(x, y);
    
    % Devuelve una función handle que evalúa el spline en un punto x
    splineFunc = @(xx) evalCubicSpline(xx, x, a, b, c, d);
end

function val = evalCubicSpline(xx, xs, as, bs, cs, ds)
    % Inicializa el valor de salida
    val = zeros(size(xx));
    
    % Para cada punto xx, encuentra el subintervalo correspondiente y evalúa el spline
    for i = 1:length(xx)
        % Encuentra el intervalo j tal que x_j <= xx < x_{j+1}
        j = find(xs(1:end-1) <= xx(i) & xs(2:end) >= xx(i), 1, 'last');
        if isempty(j)  % Si xx(i) está fuera del rango de xs
            if xx(i) < xs(1)
                j = 1;
            elseif xx(i) > xs(end)
                j = length(xs) - 1;
            end
        end
        
        % Asegúrate de que j esté dentro de los límites válidos para los vectores de coeficientes
        j = min(j, length(as));
        
        % Calcula la diferencia x - x_j
        dx = xx(i) - xs(j);
        
        % Evalúa el polinomio cúbico para el intervalo j
        val(i) = as(j) + bs(j)*dx + cs(j)*dx.^2 + ds(j)*dx.^3;
    end
end
