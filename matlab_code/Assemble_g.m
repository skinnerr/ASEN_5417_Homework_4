function [diag, sub, sup, rhs] = Assemble_g( h, BC, F, g, th )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    N = length(F);
    
    diag_range = 2:N-1;
     sub_range = 3:N-1;
     sup_range = 2:N-2;
    
    diag = (-2/h^2) - 2 *  g(diag_range);
     sub = ( 1/h^2) - 3 *  F(sub_range) / (2*h);
     sup = ( 1/h^2) + 3 *  F(sup_range) / (2*h);
     rhs =          - 1 * th(diag_range);
    
    % Account for boundary conditions, even though g0 = gf = 0.
    rhs(1)   = rhs(1)   - ((1/h^2) - 3 * F(diag_range(1))   / (2*h)) * BC.g0;
    rhs(end) = rhs(end) - ((1/h^2) + 3 * F(diag_range(end)) / (2*h)) * BC.gf;

end