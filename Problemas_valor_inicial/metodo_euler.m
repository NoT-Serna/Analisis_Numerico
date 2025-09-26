function [t, y] = euler(f, t0, y0, tf, h)
    % Método de Euler para resolver PVI
    % f   -> función @(t,y)
    % t0  -> tiempo inicial
    % y0  -> condición inicial
    % tf  -> tiempo final
    % h   -> paso de integración
    
    N = round((tf - t0)/h);  % número de pasos
    t = t0:h:tf;             % vector de tiempos
    y = zeros(1, N+1);       % inicializar solución
    y(1) = y0;               % condición inicial
    
    for n = 1:N
        y(n+1) = y(n) + h*f(t(n), y(n));
    end
end
