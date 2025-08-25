function yi = interm_value(x, y, x0)
    n = length(x);

    % Verificar rango
    if x0 < x(1) || x0 > x(end)
        error('El valor está fuera del rango de datos.');
    end

    % Localizar intervalo
    idx = find(x <= x0, 1, 'last');

    % Primer intervalo → polinomio cuadrático con los 3 primeros puntos
    if idx <= 2
        xi = x(1:3);
        fi = y(1:3);

    % Último intervalo → polinomio cuadrático con los 3 últimos puntos
    elseif idx >= n-1
        xi = x(end-2:end);
        fi = y(end-2:end);

    % Intervalos intermedios → polinomio cúbico (4 puntos)
    else
        xi = x(idx-1:idx+2);
        fi = y(idx-1:idx+2);
    end

    % Evaluar polinomio de Lagrange
    yi = lagrange_interp(xi, fi, x0);
end

function yi = lagrange_interp(x, y, x0)
    m = length(x);
    yi = 0;
    for i = 1:m
        L = 1;
        for j = 1:m
            if j ~= i
                L = L * (x0 - x(j)) / (x(i) - x(j));
            end
        end
        yi = yi + y(i) * L;
    end
end
%{
------------USO-----------------
Ej: x= 1...10 f(x) = ln(x)
x = 1:10;
y = log(x);

x0 = 6.4;
yi = interm_value(x, y, x0);
y_real = log(x0);

Ea = abs(y_real - yi);               % error absoluto

den = abs(y_real);
if den == 0
    Er = NaN;                        % no definido (f(x0)=0)
    % Alternativa: Er = Ea / max(den, 1e-12);
    % Alternativa: Er = Ea / max(abs(yi), 1e-12);  % relativo al aproximado
else
    Er = Ea / den;                   % error relativo
end

fprintf('x0 = %.2f\n', x0);
fprintf('Valor real        = %.6f\n', y_real);
fprintf('Valor interpolado = %.6f\n', yi);
fprintf('Error absoluto    = %.6e\n', Ea);
fprintf('Error relativo    = %.6e (%.6f %%)\n', Er, 100*Er);

%}
