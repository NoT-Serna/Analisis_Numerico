%Se han recolectado datos para una sección transversal de un río donde
% y = distancia desde la orilla
% H = profundidad del río
% U = velocidad del río

% Note que:
% Ac = sección transversal = integral de 0 a y H(y) dy
% Q = Caudal del río = integral de 0 a y de H(y)*U(y) dy

% Datos:
y = [0,1,3,5,7,8,10];
H = [0,1.5,3,3.5,3.2,2,0];
U = [0,0.1,0.12,0.25,0.3,0.15,0];

% Utilizando regla del trapecio compuesta
% El área de la sección transversal del río Ac:
spline = createCubicSpline(y,H);
I = trap_compuesto(spline,0,10,1000);
disp(['El área de la sección transversal del río : ', num2str(I)]);

% La profundidad promedio:
prom = I / (y(end) - y(1));
disp(['La profundidad promedio del río : ', num2str(prom)]);

%La velocidad promedio del agua:
spline2 = createCubicSpline(y,U);
I_U = trap_compuesto(spline2, 0, 10, 1000);
vel_prom = I_U / (y(end) - y(1));
disp(['La velocidad promedio del agua : ', num2str(vel_prom)]);

% Caudal del río:
Q = H .* U;
spline3 = createCubicSpline(y,Q);
I_Q = trap_compuesto(spline3, 0, 10, 1000);
disp(['El caudal del río : ', num2str(I_Q)]);
