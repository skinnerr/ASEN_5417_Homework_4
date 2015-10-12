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
     F = nan(nPr, N);
     g = nan(nPr, N);
    th = nan(nPr, N);
    
    % Initial guesses for solution variables, based on linear approximation to Homework 3.
     F(1,:) = linspace(0, 0.50, N);
     F(2,:) = linspace(0, 0.25, N);
     g(1,:) = sin(linspace(0, 1, N)*pi) / 2;
     g(2,:) = sin(linspace(0, 1, N)*pi) / 2;
    th(1,:) = linspace(1, 0, N);
    th(2,:) = linspace(1, 0, N);
    
    % Convergence criterion.
    epsilon = 0.1;
    
    hf = figure();
    
    %%%
    % Solve our three equations (for F, g, and theta) iteratively.
    %%%
    
    % Loop through Prandtl numbers.
%     for iPr = 1:nPr
    for iPr = 1:1
        
        fprintf('Working on Pr = %i\n', Pr(iPr));
        
        % Loop until convergence criterion is met.
        norm = inf;
        iteration = 1;
%         while norm > epsilon
        while iteration < 15
        
            % Containers for previous iterations' values to determine convergence.
             F_prev =  F(iPr,:);
             g_prev =  g(iPr,:);
            th_prev = th(iPr,:);
            
            figure(hf);
            hold on;
            plot(eta, F_prev);
    
            % STEP 1: Solve the g-equation.
            
%             [a, b, c, rhs] = Assemble_g( h, BC, F_prev, g_prev, th_prev );
            [a, b, c, rhs] = Assemble_g( h, BC, F(iPr,:), g(iPr,:), th(iPr,:));
            sol = Thomas(a, b, c, rhs);
            g(iPr,:) = [BC.g0; sol; BC.gf];

            % STEP 2: Solve the theta-equation.

%             [a, b, c, rhs] = Assemble_th( h, BC, Pr(iPr), F_prev );
            [a, b, c, rhs] = Assemble_th( h, BC, Pr(iPr), F(iPr,:) );
            sol = Thomas(a, b, c, rhs);
            th(iPr,:) = [BC.th0; sol; BC.thf];

            % STEP 3: Integrate g to obtain F using the Euler method.

%             F(iPr,:) = Euler(g_prev, BC.g0, h);
            F(iPr,:) = Euler(g(iPr,:), BC.g0, h);
            
            % STEP 4: Assess convergence.

             F_norm = sum(abs( F_prev -  F(iPr,:))) / N;
             g_norm = sum(abs( g_prev -  g(iPr,:))) / N;
            th_norm = sum(abs(th_prev - th(iPr,:))) / N;

            norm = F_norm + g_norm + th_norm;

            fprintf('Iteration: %02i,  norm: %8.2e, F: %8.2e, g: %8.2e, th: %8.2e, \n', ...
                           iteration,         norm,   F_norm,   g_norm,   th_norm);
            
            iteration = iteration + 1;
        
        end
            
        figure();
        hold on;
        plot(eta, F(iPr,:), 'DisplayName', 'F');
        plot(eta, g(iPr,:), 'DisplayName', 'g');
        plot(eta,th(iPr,:), 'DisplayName', 'theta');
        title(sprintf('Iteration: %i',iteration));
        xlabel('eta');
        hleg = legend('show');
        set(hleg, 'Location', 'best');
        
    end
    
    %%%
    % Process results.
    %%%
    
    % TODO TODO TODO
    
end

function [a, b, c, rhs] = Assemble_g( h, BC, F, g, th )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %     a -- diagonal
    %     b -- sub-diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
    N = length(F);
    
    a_range = 2:N-1;
    b_range = 2:N-2;
    c_range = 3:N-1;
    
    a = (-2/h^2) - 2 *  g(a_range);
    b = ( 1/h^2) - 3 *  F(b_range) / (2*h);
    c = ( 1/h^2) + 3 *  F(c_range) / (2*h);
    rhs =        - 1 * th(a_range);
    
    % Account for boundary conditions, even though g0 = gf = 0.
    rhs(1)   = rhs(1)   - ((1/h^2) - 3 * F(a_range(1))   / (2*h)) * BC.g0;
    rhs(end) = rhs(end) - ((1/h^2) + 3 * F(a_range(end)) / (2*h)) * BC.gf;

end

function [a, b, c, rhs] = Assemble_th( h, BC, Pr, F )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the theta-system.
    %     a -- diagonal
    %     b -- sub-diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %%%
    
    N = length(F);
    
    a_range = 2:N-1;
    b_range = 2:N-2;
    c_range = 3:N-1;
    
    a = (-2/h^2) * ones(1, length(a_range));
    b = ( 1/h^2) - 3 * Pr * F(b_range) / (2*h);
    c = ( 1/h^2) + 3 * Pr * F(c_range) / (2*h);
    rhs = zeros(length(a_range),1);
    
    % Account for boundary conditions, even though thf = 0.
    rhs(1)   = rhs(1)   - ((1/h^2) - 3 * Pr * F(a_range(1))   / (2*h)) * BC.th0;
    rhs(end) = rhs(end) - ((1/h^2) + 3 * Pr * F(a_range(end)) / (2*h)) * BC.thf;
    
end

function [ F ] = Euler( g, F0, delta )
    
    %%%%%%
    % Integrates the data set g using the Euler method, given the initial value F0.
    %%%
    
    F = nan(length(g),1);
    F(1) = F0;
    for i = 1:(length(g)-1)
        F(i+1) = F(i) + delta * g(i);
    end
    
end






