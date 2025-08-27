function p = biseccion_array(x_data, y_data, a, b, TOL, N0)
    % Verificar que x_data esté ordenado
    if any(diff(x_data) <= 0)
        error('x_data debe estar en orden ascendente sin duplicados.')
    end

    % Verificar que a y b estén dentro del rango de x_data
    if a < x_data(1) || b > x_data(end)
        error('a y b deben estar dentro del rango de x_data.')
    end

    i = 1;
    FA = interp1(x_data, y_data, a);
    p = a; % Inicializar p

    while i <= N0
        p_old = p;
        p = a + (b - a)/2;
        FP = interp1(x_data, y_data, p);

        if p ~= 0
            Er = abs((p - p_old)/p);
        else
            Er = abs(b - a);
        end

        if FP == 0 || Er < TOL
            fprintf('La raíz aproximada es: %.5f\n', p)
            fprintf('Número de iteraciones: %d\n', i)
            return
        end

        i = i + 1;

        if FA * FP > 0
            a = p;
            FA = FP;
        else
            b = p;
        end
    end

    error('El método falló después del número máximo de iteraciones')
end


&}
----------------USO----------------------
x = linspace(0, 10, 100);
y = sin(x); % f(x) = sin(x)

a = 3;  % Intervalo donde sin(x) cambia de signo
b = 4;
TOL = 1e-6;
N0 = 50;

p = biseccion_array(x, y, a, b, TOL, N0);


&{
