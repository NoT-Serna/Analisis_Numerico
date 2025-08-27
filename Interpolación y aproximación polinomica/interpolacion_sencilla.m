function y = my_interp(x_data, y_data, x)
    % Linear interpolation (custom)
    if x < x_data(1) || x > x_data(end)
        error('x is out of bounds')
    end

    % Find the interval [x_i, x_{i+1}] such that x_i <= x <= x_{i+1}
    for i = 1:length(x_data) - 1
        if x >= x_data(i) && x <= x_data(i + 1)
            % Linear interpolation formula
            t = (x - x_data(i)) / (x_data(i+1) - x_data(i));
            y = y_data(i) + t * (y_data(i+1) - y_data(i));
            return
        end
    end

    error('Interpolation failed')
end
