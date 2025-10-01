%Modelo de proporcionalidad donde tenemos
b = 0.02;
r = 0.1;
d = 0.015;
% Condiciones iniciales
t0 = 0;
p0 = 0.01;
tf = 50;
h = 1;
tspan = t0:h:tf;

%Ecuacion diferencial
dPdt = @(t,p) r*b*(1-p);

%Resolver con Runge-Kutta 4
[t_rk4,p_rk4] = runge_kutta_4(dPdt,t0,p0,h,tf);

%Grafica
plot(t_rk4, p_rk4);
xlabel('Time');
ylabel('Population');
title('Population Growth using Runge-Kutta 4 Method');
grid on;

%Podemos evidenciar un crecimiento lineal donde en t = 50 la poblacion
%disminuye un 0.1
