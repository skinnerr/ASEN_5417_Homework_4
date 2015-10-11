function [] = Problem_1()

    %%%%%%
    % Solves the differential equations governing free convection along a vertical plate
    % using central differences for the two coupled equations and the RK4 method for the
    % uncoupled equation.
    %
    % Ryan Skinner, October 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Solution domain.
    eta0 = 0;
    etaf = 10;
    N = 101;
    eta = linspace(eta0, etaf, N)';
    h = eta(2) - eta(1);
    
    % Prandtl number.
    Pr = [1, 10];
    nPr = length(Pr);
    
    % Boundary conditions (0 = initial, f = final).
    BC.F0  = 0;
    BC.g0  = 0;
    BC.gf  = 0;
    BC.th0 = 1;
    BC.thf = 0;
    
    % Solution variables (indexed by Prandtl number and then time step).
     F = zeros(nPr, N);
     g = zeros(nPr, N);
    th = zeros(nPr, N);
    
    % Initial guesses for solution variables, based on linear approximation to Homework 3.
     F(1,:) = linspace(0, 0.50, N);
     F(2,:) = linspace(0, 0.25, N);
     g(1,:) = linspace(0,    0, N);
     g(2,:) = linspace(0,    0, N);
    th(1,:) = linspace(1,    0, N);
    th(2,:) = linspace(1,    0, N);
    
    % Convergence criterion.
    epsilon = 0.1;
    
    %%%
    % Solve our three equations (for F, g, and theta) iteratively.
    %%%
    
    % Loop through Prandtl numbers.
    for iPr = 1:nPr
        
        % Loop until convergence criterion is met.
        norm = inf;
        iteration = 0;
        while norm < epsilon
            
            iteration = iteration + 1;
            
            % Containers for previous iterations' values to determine convergence.
             F_prev =  F(iPr,:);
             g_prev =  g(iPr,:);
            th_prev = th(iPr,:);
    
            % STEP 1: Solve the g-equation.
            
            [a, b, c, rhs] = Assemble_g( h, BC, F(iPr,:), g(iPr,:), th(iPr,:) );
            g(iPr,:) = Thomas(a, b, c, rhs);

            % STEP 2: Solve the theta-equation.

            [a, b, c, rhs] = Assemble_th( h, BC, Pr(iPr), F(iPr,:) );
            th(iPr,:) = Thomas(a, b, c, rhs);

            % STEP 3: Integrate g to obtain F using the Euler method.

            F(iPr,:) = Euler(g(iPr,:), BC.g0, h);
            
            % STEP 4: Assess convergence.

             F_norm = sum( F_prev -  F(iPR,:)) / N;
             g_norm = sum( g_prev -  g(iPR,:)) / N;
            th_norm = sum(th_prev - th(iPR,:)) / N;

            norm = F_norm + g_norm + th_norm;

            fprintf('Iteration: %02i,  norm: %10.2e', iteration, norm);
        
        end
        
    end
    
end

function [a, b, c, rhs] = Assemble_g( h, BC, F, g, th )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %     a -- sub-diagonal
    %     b -- diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
    N = length(F);
    
    a_range = 2:N-2;
    b_range = 2:N-1;
    c_range = 3:N-1;
    
    a = ( 1/h^2) - 3 *  F(a_range) / (2*h);
    b = (-2/h^2) - 2 *  g(b_range);
    c = ( 1/h^2) + 3 *  F(c_range) / (2*h);
    rhs =        - 1 * th(b_range);
    
    % Account for boundary conditions, even though g0 = gf = 0.
    rhs(1)   = rhs(1)   - ((1/h^2) - 3 * F(b_range(1))   / (2*h)) * BC.g0;
    rhs(end) = rhs(end) - ((1/h^2) + 3 * F(b_range(end)) / (2*h)) * BC.gf;

end

function [a, b, c, rhs] = Assemble_th( h, BC, Pr, F )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the theta-system.
    %     a -- sub-diagonal
    %     b -- diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
    N = length(F);
    
    a_range = 2:N-2;
    b_range = 2:N-1;
    c_range = 3:N-1;
    
    a =  1      + 3 * Pr * F(a_range) / (2*h);
    b = -2;
    c = (1/h^2) - 3 * Pr * F(c_range) / (2*h);
    rhs = zeros(length(b_range),1);
    
    % Account for boundary conditions, even though thf = 0.
    rhs(1)   = rhs(1)   - ( 1      + 3 * Pr * F(b_range(1))   / (2*h)) * BC.th0;
    rhs(end) = rhs(end) - ((1/h^2) - 3 * Pr * F(b_range(end)) / (2*h)) * BC.thf;
    
end

function [ F ] = Euler( g, F0, delta )
    
    %%%%%%
    % Integrates the data set g using the Euler method, given the initial value F0.
    %%%
    
    F = zeros(length(g),1);
    F(1) = F0;
    for i = 1:(length(g)-1)
        F(i+1) = F(i) + delta * g(i);
    end
    
end