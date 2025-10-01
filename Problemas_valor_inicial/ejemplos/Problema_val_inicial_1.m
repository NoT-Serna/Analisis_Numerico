%Considerando un circuito tenemos que ver el cambio de la corriente (i) en
%un tiempo (t)

%Parámetros del circuito son:
C = 0.3;
R = 1.4;
L = 1.7;
% Fuente de voltaje
e = @(t) exp(-0.067*t)* sin(2*t-pi);
%Necesitamos la primera y segunda derivada de e
ePrime = @(t) 0.067*exp(-0.067*t)*sin(2*t)-2*exp(-0.067*t)*cos(2*t);
ePrime2 = @(t) 3.995511*exp(-0.067*t)*sin(2*t)+ 0.268*exp(-0.067*t)*cos(2*t);

%Condiciones iniciales
t0 = 0;
i0 = 0;
tf = 10;
h = 0.1;
tspan = t0:h:tf;

%Ecuación diferencial
dIdt = @(t, i) C*ePrime2(t)+(1/R)*ePrime(t)+(1/L)*e(t);

%Resolver con Runge-kutta 4
[t_rk4, i_rk4] = runge_kutta_4(dIdt,t0,i0,h,tf);
% Graficar la corriente obtenida
figure;
plot(t_rk4, i_rk4);
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
title('Cambio de Corriente en el Circuito');
grid on;
