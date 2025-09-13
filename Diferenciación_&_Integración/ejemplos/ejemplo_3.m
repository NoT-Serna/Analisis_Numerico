%Se desea determinar que tan lejos ha caído un paracaidista despúes de un
%tiempo t = 9.9334

% los datos que tenemos son:
t = 9.9334; % tiempo en segundos
g = 9.8; % aceleración debida a la gravedad en m/s^2
m = 68.1; % masa en kilogramos
c = 12.5; % constante en kg/s

% luego la formulade distancia esta dada por:
f = @(t) (g*m)/c * (1-exp((-c/m).*t));
% Sabemos que analiticamente la función es d = 289.43515 ahora veamos
% nuestra arpoximación
d = 2289.43515;

% Nuestro resultado usando la regla del trapecio compuesto:
I = trap_compuesto(f,0,t,1000);
%Aproximacion de 5 decimales para I
res = sprintf('%.5f', I);


disp(['La distancia del paracaidista es de : ', num2str(I)]);


% Diferencia entre d y nuestra aproximacion
error = abs(d - I);
fprintf('El error entre la distancia analítica y la aproximación es: %.5f\n', error);
% La distancia del paracaidista es de : 286.45 m
