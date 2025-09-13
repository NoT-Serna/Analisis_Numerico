%La masa total de una barra con densidad variable está dada por 
% m = integral de 0 a L de p(x) * Ac(x) dx
% Donde 
% L = Longitud de la barra en cm
% Ac(x) = Área de la sección transversal en cm^2 (constante en este caso)
% p = Densidad variable en g/cm^3

%Datos
L = 10;
x = [0,2,3,4,6,8,10];
p = [4.00,3.95,3.89,3.80,3.60,3.41,3.30];
Ac = [100,103,106,110,120,133,150];

%Generar funciones (splines) para operar en la integral
spline = createCubicSpline(x,p);
spline2 = createCubicSpline(x,Ac);
f = @(x) spline(x).*spline2(x);
% Masa total:
m = trap_compuesto(f,0,L,1000);
% Mostrar el resultado de la masa total
disp(['La masa total de la barra es: ', num2str(m), ' g']);

%Conversion a kg
p_kg = p * 1000;
Ac_kg = Ac * 1e-4;
spline3 = createCubicSpline(x,p_kg);
spline4 = createCubicSpline(x,Ac_kg);
f2 = @(x) spline3(x).*spline4(x);
% Masa total en kg:
m_kg = trap_compuesto(f2, 0, L, 1000);
% Mostrar el resultado de la masa total en kg
disp(['La masa total de la barra en kg es: ', num2str(m_kg), ' kg']);

