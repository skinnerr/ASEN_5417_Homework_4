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
    
    % Containers for previous iterations' values to determine convergence.
     Fprev = 0;
     gprev = 0;
    thprev = 0;
    
    % Loop through Prandtl numbers.
    for iPr = 1:nPr
        
        % Loop until convergence criterion is met.
        norm = inf;
        iteration = 0;
        while norm < epsilon
            
            iteration = iteration + 1;
             Fprev =  F(iPr,:);
             gprev =  g(iPr,:);
            thprev = th(iPr,:);
    
            % STEP 1: Solve the g-equation.
            
            [sub, diag, sup, rhs] = Assemble_g( h, Pr(iPr), F(iPr,:), g(iPr,:), th(iPr,:) );
            g(iPr,:) = Thomas(sub, diag, sup, rhs);

            % STEP 2: Solve the theta-equation.

            [sub, diag, sup, rhs] = Assemble_th( h, Pr(iPr), F(iPr,:) );
            th(iPr,:) = Thomas(sub, diag, sup, rhs);

            % STEP 3: Integrate g to obtain F using the Euler method.

            F(iPr,:) = Euler(g(iPr,:), BC.g0, h);
            
            % STEP 4: Assess convergence.

             Fnorm = sum( Fprev -  F(iPR,:)) / N;
             gnorm = sum( gprev -  g(iPR,:)) / N;
            thnorm = sum(thprev - th(iPR,:)) / N;

            norm = Fnorm + gnorm + thnorm;

            fprintf('Iteration: %02i,  norm: %10.2e', iteration, norm);
        
        end
        
    end
    
end

function [a, b, c, rhs] = Assemble_g( h, Pr, F, g, th )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %     a -- sub-diagonal
    %     b -- diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
    N = length(F);
    a = zeros(N-2);
    b = zeros(N-2);
    c = zeros(N-2);

end

function [a, b, c, rhs] = Assemble_th( h, Pr, F )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the theta-system.
    %     a -- sub-diagonal
    %     b -- diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
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