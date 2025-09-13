%Dimensiones del lago
x = [0,3,6,9,12,15,18,21,24,27,30,33,36,39,42];
y = [0,6,7,8,10,9,7,4,6,4,5,4,2,1,0];
% Aproixmar a una función con splines cúbicos
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

% Como queremos ver el volumen de agua necesario es decir s^2 *pi evaluado
% en Simpson
% Calcular el volumen de agua usando la fórmula del volumen
f = @(t) (spline(t)).^2;
Volumen = pi * Simpson(0,42,1000000,f);
disp(['El volumen de agua necesario es: ', num2str(Volumen)]);
% Se necesitan 4631.0579 metros cúbicos de agua para rellenarlo