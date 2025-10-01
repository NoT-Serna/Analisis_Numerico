% Parámetros
m = 0.11;
g = 9.8;
k = 0.002;

% Fuerza de resistencia
Fr = @(v) -k*v*abs(v);

% Ecuación diferencial
f = @(t,v) (-m*g + Fr(v)) / m;

% Condiciones iniciales
t0 = 0;
n = 10;     % hasta 1 segundo
h = 0.1;    % paso
v0 = 8;

% Llamada a Runge-Kutta 4
[t,v] = runge_kutta_4(f,t0,v0,h,n);

% Resultados
% Resultados (solo desde 0.1 hasta 1.0)
disp('  t (s)     v (m/s)')
disp('---------------------')
for i = 2:length(t)   % empieza en 2 porque t(1)=0
    fprintf('%6.2f   %10.6f\n', t(i), v(i));
end


% Gráfica
plot(t, v, 'o-','LineWidth',1.2)
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
title('Velocidad usando Runge-Kutta 4')
grid on
