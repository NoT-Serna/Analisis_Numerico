function I = trap_compuesto(f, a, b, n)
% I = trap_compuesto(f,a,b,n) calcula la integral aproximada de f en [a,b]
% usando la regla del trapecio compuesto con n subintervalos.
h = (b - a) / n;
x = a : h : b;
y = f(x);
I = h * (0.5 * y(1) + sum(y(2:end-1)) + 0.5 * y(end));
end