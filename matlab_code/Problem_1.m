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
    % Calculate the solutions from Homework 3, for Pr = {1, 10}.
    %%%
    
    [T1, Y1, T10, Y10] = Problem_1_Shooting(false);
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Solution domain.
    eta0 = 0;
    etaf = 10;
    N = 501;
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
    epsilon = 1e-5;
    
    %%%
    % Solve our three equations (for F, g, and theta) iteratively.
    %%%
    
    % Loop through Prandtl numbers.
    for iPr = 1:nPr
        
        fprintf('Working on Pr = %i\n', Pr(iPr));
        
        % Loop until convergence criterion is met.
        norm = inf;
        iteration = 1;
        while norm > epsilon
        
            % Containers for previous iterations' values to determine convergence.
             F_prev =  F(iPr,:);
             g_prev =  g(iPr,:);
            th_prev = th(iPr,:);
    
            %%%
            % STEP 1: Solve the g-equation.
            %%%
            
            [diag, sub, sup, rhs] = Assemble_g( h, BC, F_prev, g_prev, th_prev );
            sol = Thomas(diag, sub, sup, rhs);
            g(iPr,:) = [BC.g0; sol; BC.gf];

            %%%
            % STEP 2: Solve the theta-equation.
            %%%

            [diag, sub, sup, rhs] = Assemble_th( h, BC, Pr(iPr), F_prev );
            sol = Thomas(diag, sub, sup, rhs);
            th(iPr,:) = [BC.th0; sol; BC.thf];

            %%%
            % STEP 3: Integrate g to obtain F using the Euler method.
            %%%

            F(iPr,:) = Euler(g_prev, BC.g0, h);
            
            %%%
            % STEP 4: Assess convergence.
            %%%

             F_norm = sum(abs( F_prev -  F(iPr,:))) / N;
             g_norm = sum(abs( g_prev -  g(iPr,:))) / N;
            th_norm = sum(abs(th_prev - th(iPr,:))) / N;

            norm = F_norm + g_norm + th_norm;

            fprintf('Iteration: %02i,  norm: %8.2e, F: %8.2e, g: %8.2e, th: %8.2e\n', ...
                           iteration,         norm,   F_norm,   g_norm,   th_norm);
            
            iteration = iteration + 1;
        
        end
        
    end
    
    %%%
    % Process results.
    %%%
    
    hF  = figure();
    hg  = figure();
    hth = figure();
    
    for iPr = 1:length(Pr)

        if Pr(iPr) == 1
            T = T1;
            Y = Y1;
        elseif Pr(iPr) == 10
            T = T10;
            Y = Y10;
        else
            error('Prandtl number %.2f not supported.',Pr(iPr));
        end

        % Plot F

        figure(hF);
        hold on;
        plot(eta, F(iPr,:),     'DisplayName',sprintf('Pr = %i (2CD)',Pr(iPr)));
        plot(  T,   Y(:,1),'--','DisplayName',sprintf('Pr = %i (HW3)',Pr(iPr)));
        xlabel('\eta');
        ylabel('F');
        
        % Plot F'
        
        figure(hg);
        hold on;
        plot(eta, g(iPr,:),     'DisplayName',sprintf('Pr = %i (2CD)',Pr(iPr)));
        plot(  T,   Y(:,2),'--','DisplayName',sprintf('Pr = %i (HW3)',Pr(iPr)));
        xlabel('\eta');
        ylabel('F''');
        
        % Plot theta
        
        figure(hth);
        hold on;
        plot(eta,th(iPr,:),     'DisplayName',sprintf('Pr = %i (2CD) ',Pr(iPr)));
        plot(  T,   Y(:,4),'--','DisplayName',sprintf('Pr = %i (HW3)',Pr(iPr)));
        xlabel('\eta');
        ylabel('\theta');

        % Determine F''(0) and theta'(0) using second-order forward differences.
        Fpp = (-3* g(iPr,1) + 4* g(iPr,2) -    g(iPr,3)) / (2*h);
        thp = (-3*th(iPr,1) + 4*th(iPr,2) -   th(iPr,3)) / (2*h);
        
        fprintf('Pr = %2i: F''''(0) = %7.4f, th''(0) = %7.4f\n', Pr(iPr), Fpp, thp);

    end
    
    for h = [hF, hg, hth]
        figure(h);
        hleg = legend('show');
        set(hleg, 'Location', 'eastoutside');
    end
    
end










