function d = derivada3_NR(y, h)
    % y: vector de datos igualmente espaciados
    % h: paso
    n = length(y);
    d = zeros(1,n);
    
    % Adelante en el primer punto
    d(1) = (-3*y(1) + 4*y(2) - y(3)) / (2*h);
    
    % Centrada en los puntos interiores
    for i = 2:n-1
        d(i) = (y(i+1) - y(i-1)) / (2*h);
    end
    
    % Atrás en el último punto
    d(n) = (3*y(n) - 4*y(n-1) + y(n-2)) / (2*h);
end

function d = derivada5_NR(y, h)
    % y: vector de datos igualmente espaciados
    % h: paso
    n = length(y);
    d = zeros(1,n);
    
    % Fórmula hacia adelante en los primeros dos puntos
    d(1) = (-25*y(1) + 48*y(2) - 36*y(3) + 16*y(4) - 3*y(5)) / (12*h);
    d(2) = (-25*y(2) + 48*y(3) - 36*y(4) + 16*y(5) - 3*y(6)) / (12*h);
    
    % Fórmula centrada en los puntos interiores
    for i = 3:n-2
        d(i) = (y(i-2) - 8*y(i-1) + 8*y(i+1) - y(i+2)) / (12*h);
    end
    
    % Fórmula hacia atrás en los dos últimos puntos
    d(n-1) = (25*y(n-1) - 48*y(n-2) + 36*y(n-3) - 16*y(n-4) + 3*y(n-5)) / (12*h);
    d(n)   = (25*y(n)   - 48*y(n-1) + 36*y(n-2) - 16*y(n-3) + 3*y(n-4)) / (12*h);
end


%{
---------USO---------------

x = linspace(0,2*pi,101);
h = x(2)-x(1);
y = sin(x);

d3 = derivada3_NR(y,h);
d5 = derivada5_NR(y,h);

plot(x,cos(x),'k-',x,d3,'r--',x,d5,'b-.')
legend('Exacta','3 puntos','5 puntos')

%}
