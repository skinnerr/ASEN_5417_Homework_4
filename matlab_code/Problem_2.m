function [] = Problem_2()

    %%%%%%
    % Solves the differential equation (y'' + y/4 = 0) using second-order central
    % differences with unequally-spaced grid points and LU-decomposition.
    %
    % Ryan Skinner, October 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Solution domain: the closed interval [-1, 1].
    N = 51;
    th = flip(cos(pi*(0:N-1)/(N-1)));
    
    % Boundary conditions (0 = initial, f = final).
    BC.y0 = 0;
    BC.yf = 2;
    
    % Assemble central difference matrix equation for y.
    [a, b, c, rhs] = Assemble_y(th, BC);
    
    % Solve matrix equation.
%     y = LUDecomp(a, b, c, rhs);
    y = Thomas(a, b, c, rhs);
    y = [BC.y0; y; BC.yf];
	
    % Compute analytical solution.
    y_exact = (2/sin(1)) * sin((th+1)/2)';
    
    figure();
    hold on;
    plot(th, y_exact, 'LineStyle', '-', 'DisplayName', 'Analytical Solution');
    plot(th,       y, 'o', 'DisplayName', '2nd-order Central Differences');
    xlabel('\theta');
    hleg = legend('show');
    set(hleg, 'Location', 'southeast');
    
    figure();
    hold on;
    plot(th, y-y_exact, '-o');
    xlabel('\theta');
    ylabel('Point-Wise Error');
    
end

function [a, b, c, rhs] = Assemble_y( th, BC )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the y-system.
    %     a -- diagonal
    %     b -- sub-diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
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
    
    a = 1/4 + D;
    b = C(2:end);
    c = B(1:end-1);
    rhs = zeros(N-2,1);
    
    % Account for boundary conditions.
    rhs(1)   = rhs(1)   - C(1)   * BC.y0;
    rhs(end) = rhs(end) - B(end) * BC.yf;

end















