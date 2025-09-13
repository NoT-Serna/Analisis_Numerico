% Determinar los arrays para el cálculo del área
x = [0,100,200,300,400,500,600,800,900,1000];
y = [125,125,120,112,90,95,88,75,35,0];
% Aproximar los datos a una función polinomica con spline cúbicos:
spline = createCubicSpline(x,y);
xx = linspace(min(x),max(x),1000);
yy = spline(xx); % Evalua el spline en xx
%Ver la gráfica sobre la vemos cual es el área a integrar:
plot(xx, yy, 'b-', 'LineWidth', 1.5);   % spline en azul
hold on;
plot(x, y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % puntos originales en rojo
legend('Spline cúbico', 'Puntos iniciales');
grid on;
hold off;

%Usar método de Simpson para integrar
s = Simpson(0,1000,1000,spline);
% Mostrar el resultado de la integración
disp(['El área bajo la curva es: ', num2str(s)]);




