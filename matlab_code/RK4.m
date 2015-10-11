function [ T, Y ] = RK4( odefun, tspan, N, y0 )

    %%%%%%
    % Solves a differential equation using the RK4 method.
    % INPUTS: odefun -- function handle to the system of odes (as in ode45)
    %          tspan -- vector specifying the interval of integration
    %              N -- number of time steps
    %             t0 -- vector of initial conditions
    %
    % Ryan Skinner, September 2015
    %%%
    
    % Number of equations to integrate.
    neq = length(y0);
    
    % Number of preliminary values calculated by RK4 method.
            nk = 4;
       t_coeff = [0, 1, 1, 1]' / 2;
       k_coeff = [0, 1, 1, 2]' / 2;
    ynp1_coeff = [1, 2, 2, 1]' / 6;

    % Time-like variable at which to evaluate solution, and step size.
    T = linspace(tspan(1), tspan(2), N)';
    h = T(2) - T(1);
    
    % Solution vector with initial conditions.
    Y = zeros(N, neq);
    Y(1,:) = y0;
    
    %%%
    % Perform integration using RK4.
    %%%
    
    % Loop over time steps.
    for n = 1:(length(T)-1)
        k = zeros(neq, nk+1);
        % Loop over Runge-Kutta intermediary values k_1, k_2, ...
        for i = (1:nk)+1
            dy = odefun( T(n)    + t_coeff(i-1) * h,       ...
                         Y(n,:)' + k_coeff(i-1) * k(:,i-1) );
            k(:,i) = h * dy; 
        end
        Y(n+1,:) = Y(n,:) + (k(:,2:nk+1) * ynp1_coeff)';
    end
    
end