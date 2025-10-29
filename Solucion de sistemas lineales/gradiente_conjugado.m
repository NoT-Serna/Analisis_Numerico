function [x_sol, vect_residual] = gradiente_conjugado(A, b, C, n, tol)
% Método del Gradiente Conjugado Preacondicionado
% -------------------------------------------------
% Resuelve el sistema A*x = b donde A es simétrica y definida positiva
% utilizando una matriz de preacondicionamiento C (por ejemplo, diag(A)).
%
% Entradas:
%   A   -> matriz cuadrada (simétrica y definida positiva)
%   b   -> vector del lado derecho
%   C   -> matriz de preacondicionamiento (idealmente SPD)
%   n   -> número máximo de iteraciones
%   tol -> tolerancia de convergencia
%
% Salidas:
%   x_sol         -> solución aproximada
%   vect_residual -> norma del residuo en cada iteración

    % Inicialización
    x = zeros(size(b));          % Aproximación inicial x0 = 0
    r = b - A * x;               % Residuo inicial
    z = C \ r;                   % Resolver C*z = r  (en lugar de inv(C)*r)
    p = z;
    vect_residual = zeros(n,1);
    vect_residual(1) = norm(r);

    for k = 1:n
        Ap = A * p;
        alpha = (r' * z) / (p' * Ap);
        x = x + alpha * p;
        r_new = r - alpha * Ap;

        vect_residual(k) = norm(r_new);
        if vect_residual(k) < tol
            break;
        end

        z_new = C \ r_new;                     % z = C^{-1} * r_new
        beta = (r_new' * z_new) / (r' * z);
        p = z_new + beta * p;

        % Actualización
        r = r_new;
        z = z_new;
    end

    % Salidas finales
    x_sol = x;
    vect_residual = vect_residual(1:k);
end


%{
%-------USO-----------%
A = [4 1; 1 3];
b = [1; 2];
C = diag(diag(A));   % Preacondicionador simple
n = 100;
tol = 1e-6;

[x_sol, res] = gradiente_conjugado(A, b, C, n, tol);
disp(x_sol)
plot(res, '-o')
xlabel('Iteración')
ylabel('Norma del residuo')
title('Convergencia del método del gradiente conjugado preacondicionado')
grid on

%}
