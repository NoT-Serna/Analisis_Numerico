function R = romberg(f,a,b,n)
    % f : función a integrar
    % [a,b] : intervalo
    % n : número de niveles de Romberg (ej. 4, 5, etc.)
    
    R = zeros(n,n);
    h = b - a;
    
    % Primer trapecio
    R(1,1) = 0.5 * h * (f(a) + f(b));
    
    % Romberg iterations
    for k = 2:n
        h = h/2;
        % Sumar los puntos intermedios nuevos (2^(k-2) puntos)
        subtotal = 0;
        for i = 1:2^(k-2)
            x = a + (2*i - 1)*h;
            subtotal = subtotal + f(x);
        end
        
        % Trapecio refinado
        R(k,1) = 0.5*R(k-1,1) + h*subtotal;
        
        % Extrapolación de Richardson
        for j = 2:k
            R(k,j) = R(k,j-1) + (R(k,j-1) - R(k-1,j-1))/(4^(j-1) - 1);
        end
    end
end
