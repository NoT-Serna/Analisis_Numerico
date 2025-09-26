function [t, y] = runge_kutta_4(f, t0, y0, h, n)
    % Inputs:
    % f  - function handle, e.g., @(t, y) y - t^2 + 1
    % t0 - initial time
    % y0 - initial value (can be scalar or vector)
    % h  - step size
    % n  - number of steps

    % Initialize output arrays
    t = zeros(n+1, 1);
    y = zeros(n+1, length(y0));

    % Initial conditions
    t(1) = t0;
    y(1, :) = y0;

    % Runge-Kutta 4 loop
    for i = 1:n
        k1 = f(t(i),          y(i, :)')';
        k2 = f(t(i) + h/2,    y(i, :)' + h/2 * k1')';
        k3 = f(t(i) + h/2,    y(i, :)' + h/2 * k2')';
        k4 = f(t(i) + h,      y(i, :)' + h   * k3')';

        y(i+1, :) = y(i, :) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        t(i+1) = t(i) + h;
    end
end

%{
Usage:

% Define the differential equation dy/dt = y - t^2 + 1
f = @(t, y) y - t^2 + 1;

% Initial condition
t0 = 0;
y0 = 0.5;

% Step size and number of steps
h = 0.2;
n = 10;

% Run RK4
[t, y] = runge_kutta_4(f, t0, y0, h, n);

% Plot result
plot(t, y, 'o-')
xlabel('t')
ylabel('y')
title('RK4 Solution')
grid on


%}
