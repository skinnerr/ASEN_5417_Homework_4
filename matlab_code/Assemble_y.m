function [diag, sub, sup, rhs] = Assemble_y( th, BC )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the y-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    N = length(th);
    
    B = nan(N-2,1);
    C = nan(N-2,1);
    D = nan(N-2,1);
    for i = 2:N-1
        B(i-1) =   2 / ((th(i+1) - th(i-1)) * (th(i+1) - th(i)));
        C(i-1) =   2 / ((th(i+1) - th(i-1)) * (th(i)   - th(i-1)));
        D(i-1) = (-2 / ( th(i+1) - th(i-1))) * (1/(th(i+1)-th(i)) + 1/(th(i)-th(i-1)));
    end
    
    diag = 1/4 + D;
    sub = C(2:end);
    sup = B(1:end-1);
    rhs = zeros(N-2,1);
    
    % Account for boundary conditions.
    rhs(1)   = rhs(1)   - C(1)   * BC.y0;
    rhs(end) = rhs(end) - B(end) * BC.yf;

end