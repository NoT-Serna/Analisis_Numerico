function [t, y] = rkf45(f, tspan, y0, h0, tol)
    % f: función @(t,y)
    % tspan: [t0 tf]
    % y0: condición inicial (vector fila o columna)
    % h0: paso inicial
    % tol: tolerancia de error
    
    % Coeficientes del método RKF45 (Numerical Recipes)
    a = [0       1/4     3/8     12/13   1       1/2];
    b = [...
        0           0           0           0           0
        1/4         0           0           0           0
        3/32        9/32        0           0           0
        1932/2197  -7200/2197   7296/2197   0           0
        439/216    -8           3680/513   -845/4104    0
        -8/27       2          -3544/2565   1859/4104  -11/40];
    c4 = [25/216  0   1408/2565  2197/4104  -1/5    0];
    c5 = [16/135  0   6656/12825 28561/56430 -9/50  2/55];

    t0 = tspan(1); tf = tspan(2);
    t = t0;
    y = y0(:)'; % fila

    h = h0;
    while t(end) < tf
        if t(end) + h > tf
            h = tf - t(end);
        end
        
        % Calcular k's
        k1 = f(t(end), y(end,:)')';
        k2 = f(t(end)+a(2)*h, y(end,:)+h*(b(2,1)*k1))';
        k3 = f(t(end)+a(3)*h, y(end,:)+h*(b(3,1)*k1+b(3,2)*k2))';
        k4 = f(t(end)+a(4)*h, y(end,:)+h*(b(4,1)*k1+b(4,2)*k2+b(4,3)*k3))';
        k5 = f(t(end)+a(5)*h, y(end,:)+h*(b(5,1)*k1+b(5,2)*k2+b(5,3)*k3+b(5,4)*k4))';
        k6 = f(t(end)+a(6)*h, y(end,:)+h*(b(6,1)*k1+b(6,2)*k2+b(6,3)*k3+b(6,4)*k4+b(6,5)*k5))';
        
        % Aproximación orden 4 y 5
        y4 = y(end,:) + h*(c4(1)*k1+c4(2)*k2+c4(3)*k3+c4(4)*k4+c4(5)*k5+c4(6)*k6);
        y5 = y(end,:) + h*(c5(1)*k1+c5(2)*k2+c5(3)*k3+c5(4)*k4+c5(5)*k5+c5(6)*k6);
        
        % Estimación del error
        err = norm(y5 - y4, inf);
        
        if err < tol
            % Paso aceptado
            t(end+1,1) = t(end) + h;
            y(end+1,:) = y5;
        end
        
        % Ajuste del paso
        if err == 0
            s = 2;
        else
            s = 0.84 * (tol/err)^(1/4);
        end
        h = h * min(2, max(0.1, s));
    end
end

%}
--Uso---
f = @(t,y) -2*t.*y;   % EDO: y' = -2ty
[t,y] = rkf45(f,[0 2],1,0.1,1e-6);

plot(t,y,'o-')
xlabel('t'), ylabel('y')
%{
