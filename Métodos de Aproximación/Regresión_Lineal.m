function [a0, a1, syx, r2] = regress(x, y)
    % Calcular el número de elementos
    n = length(x);

    % Inicialización de las variables de suma
    sumx = 0; 
    sumy = 0; 
    sumxy = 0; 
    sumx2 = 0; 
    st = 0; 
    sr = 0;

    % Cálculo de sumas necesarias
    for i = 1:n
        sumx = sumx + x(i);
        sumy = sumy + y(i);
        sumxy = sumxy + x(i) * y(i);
        sumx2 = sumx2 + x(i)^2;
    end

    % Cálculo de los promedios de x y y
    xm = sumx / n;
    ym = sumy / n;

    % Cálculo de los coeficientes de la recta de regresión
    a1 = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx^2);
    a0 = ym - a1 * xm;

    % Cálculo de ST y SR
    for i = 1:n
        st = st + (y(i) - ym)^2;
        sr = sr + (y(i) - a1 * x(i) - a0)^2;
    end

    % Cálculo de la desviación estándar y el coeficiente de determinación
    syx = sqrt(sr / (n - 2));
    r2 = (st - sr) / st;
end
