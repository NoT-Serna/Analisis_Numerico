function [a, b, c, d] = cubicSpline(x, a)
    n = length(x) - 1; % Número de intervalos
    h = diff(x); % Diferencias entre x consecutivos: h_i = x_{i+1} - x_i
    
    % Paso 1: Calcular los valores de alpha
    alpha = zeros(n-1, 1);
    for i = 2:n
        alpha(i-1) = (3/h(i))*(a(i+1) - a(i)) - (3/h(i-1))*(a(i) - a(i-1));
    end
    
    % Paso 2: Resolver el sistema tridiagonal para encontrar c
    l = zeros(n+1, 1);
    mu = zeros(n+1, 1);
    z = zeros(n+1, 1);
    
    l(1) = 1;
    for i = 2:n
        l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1);
        mu(i) = h(i)/l(i);
        z(i) = (alpha(i-1) - h(i-1)*z(i-1))/l(i);
    end
    
    l(n+1) = 1;
    c = zeros(n+1, 1);
    b = zeros(n, 1);
    d = zeros(n, 1);
    
    % Paso 3: Calcular b, c, d hacia atrás
    for j = n:-1:1
        c(j) = z(j) - mu(j)*c(j+1);
        b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2*c(j))/3;
        d(j) = (c(j+1) - c(j))/(3*h(j));
    end
    
    % Los coeficientes a son los mismos que los valores de a dados
    a = a(1:n); % Ajusta el tamaño de a para que coincida con b, c, d
end
