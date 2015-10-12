function [ x0, x ] = Secant1D( f, x_initials, tol )

    %%%%%%
    % Simple secant method to find root of function, f(x), based on an initial guess
    %  of the interval (x_initials). Solution is found when the relative error decreases
    %  below the tolerance, tol.
    %
    % Ryan Skinner, September 2015
    %%%
    
    % Initialize the test points.
    x = x_initials;
    y = [f(x(1)), f(x(2))];
    
    % Iterate until tolerance is met.
    relerr = inf;
    while relerr > tol
        
        % Calculate next point using the secant method.
        x(end+1) = x(end) - y(end) * (x(end) - x(end-1)) / (y(end) - y(end-1));
        
        % Calculate the next y-value.
        y(end+1) = f(x(end));
        
        % Re-calculate relative error.
        relerr = abs((x(end) - x(end-1)) / x(end-1));
        
        % Print convergence if desired.
        fprintf('x: %10.5f, err: %10.5e\n', x(end), relerr);
        if isnan(x(end))
            error('Secant method diverged.');
        end
        
    end
    
    % Return most recent guess of x.
    x0 = x(end);

end

