function [ x ] = LU_Solve( b, l, u, rhs )
    
    %%%%%%
    % Solves an LU matrix system for the solution vector x.
    %     b -- sub-diagonal of L matrix
    %     l -- diagonal of L matrix
    %     u -- super-diagonal of U matrix
    %   rhs -- right-hand side vector of matrix equation
    %     x -- solution vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    N = length(l);
    
    z = nan(N,1);
    x = nan(N,1);
    
    z(1) = rhs(1) / l(1);
    for i = 2:N
        z(i) = (rhs(i) - b(i-1) * z(i-1)) / l(i);
    end
    
    x(N) = z(N);
    for i = N-1:-1:1
        x(i) = z(i) - u(i) * x(i+1);
    end

end