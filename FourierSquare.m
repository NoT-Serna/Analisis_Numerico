function [t, f_sum] = FourierSquare(A0, T, n)
% Esta funcion grafica los primeros 'n' terminos de la serie de Fourier
% para una onda cuadrada (que inicia en +A0) y su suma total.

% 1. Crear el vector de tiempo
t = linspace(0, 4*T, 2000); 

% 2. Inicializar el vector de suma
f_sum = zeros(size(t));

% 3. Configurar la figura
figure;
hold on;
grid on;
title(['Aproximacion de Onda Cuadrada con ' num2str(n) ' Terminos']);
xlabel('Tiempo (t)');
ylabel('f(t)');

legend_strings = {}; 

% 4. Bucle para calcular y graficar cada termino
for i = 1:n
    % Indice armonico impar (k = 1, 3, 5, ...)
    k = 2*i - 1;
    
    % Calcular el termino actual
    termino_actual = (4*A0 / (k*pi)) * sin(2*pi*k*t / T);
    
    % Graficar (ya no guardamos el handle 'h')
    plot(t, termino_actual, 'r:');      
    
    % Anadir el texto de la leyenda
    legend_strings{end+1} = ['Termino k=' num2str(k)]; 
    
    % Acumular la suma
    f_sum = f_sum + termino_actual;
end

% 5. Graficar la suma total (ya no guardamos 'h_sum')
plot(t, f_sum, 'k-', 'LineWidth', 2); 
legend_strings{end+1} = 'Suma Total'; 

% 6. Ajustar la grafica
hold off;

% Llamar a legend con los textos. MATLAB los asignara
legend(legend_strings, 'Location', 'best'); 



xlim([0, 4*T]); % Limitar el eje x de 0 a 4T

end
