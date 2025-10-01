%Aproximar el tiempo requerido para que la particula reduzca si velocidad a
%v = 5m/s

%Datos
m = 10;
v0 = 10;
vf = 5;
R = @(v) -v*sqrt(v);

%Aproximaxcion de integral
t = @(u) (m./R(u));

s = Simpson(10,5,1000,t);

disp(['El tiempo requerido para que la particula disminuya a 5m/s es : ', num2str(s)]);
