function [] = Problem_2()

    %%%%%%
    % Solves the boundary value problem for free convection along a plate for Prandtl
    %  numbers of Pr = {1, 10}, using an RK4 integrator and secant-method root-finder.
    %
    % Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % Initialize root-finding functions for Prandtl numbers Pr = {1, 10}.
    rk4_Pr1  = @(x) theta10(x, 1,  1);
    rk4_Pr10 = @(x) theta10(x, 10, 1);
    o45_Pr1  = @(x) theta10(x, 1,  2);
    o45_Pr10 = @(x) theta10(x, 10, 2);
    
    % Find value of F'' using the secant method.
    fprintf('(RK4)   Pr = 1\n');
    [a.x0, a.x] = Secant1D(rk4_Pr1,  [0.60, 0.61], 1e-6);
    fprintf('(RK4)   Pr = 10\n');
    [b.x0, b.x] = Secant1D(rk4_Pr10, [0.41, 0.46], 1e-6);
    fprintf('(ode45) Pr = 1\n');
    [c.x0, c.x] = Secant1D(o45_Pr1,  [0.60, 0.61], 1e-6);
    fprintf('(ode45) Pr = 10\n');
    [d.x0, d.x] = Secant1D(o45_Pr10, [0.41, 0.46], 1e-6);
    
    % Report values of F''.
    fprintf('(RK4)   Pr = 1 :   F''''(0) = %.5f\n', a.x0);
    fprintf('(RK4)   Pr = 10:   F''''(0) = %.5f\n', b.x0);
    fprintf('(ode45) Pr = 1 :   F''''(0) = %.5f\n', c.x0);
    fprintf('(ode45) Pr = 10:   F''''(0) = %.5f\n', d.x0);
    
    % Solve differential system with optimized values of F''.
    [Ta,Ya] = RK4(  @(t, y) convection(t, y, 1),  [0,10], 500, [0,0,a.x0,1,-0.5671]);
    [Tb,Yb] = RK4(  @(t, y) convection(t, y, 10), [0,10], 500, [0,0,b.x0,1,-1.17]);
    [Tc,Yc] = ode45(@(t, y) convection(t, y, 1),  [0,10],      [0,0,c.x0,1,-0.5671]);
    [Td,Yd] = ode45(@(t, y) convection(t, y, 10), [0,10],      [0,0,d.x0,1,-1.17]);
    
    % Plot F.
    figure();
    hold on;
    plot(Ta,Ya(:,1), '-','DisplayName','(RK4) Pr = 1');
    plot(Tb,Yb(:,1), '-','DisplayName','(RK4) Pr = 10');
    plot(Tc,Yc(:,1),'--','DisplayName','(ode45) Pr = 1');
    plot(Td,Yd(:,1),'--','DisplayName','(ode45) Pr = 10');
    xlabel('\eta');
    ylabel('F');
    hleg = legend('show');
    set(hleg,'Location','northeastoutside');

    % Plot F'.
    figure();
    hold on;
    plot(Ta,Ya(:,2), '-','DisplayName','(RK4) Pr = 1');
    plot(Tb,Yb(:,2), '-','DisplayName','(RK4) Pr = 10');
    plot(Tc,Yc(:,2),'--','DisplayName','(ode45) Pr = 1');
    plot(Td,Yd(:,2),'--','DisplayName','(ode45) Pr = 10');
    xlabel('\eta');
    ylabel('F''');
    hleg = legend('show');
    set(hleg,'Location','northeastoutside');
    
    % Plot theta.
    figure();
    hold on;
    plot(Ta,Ya(:,4), '-','DisplayName','(RK4) Pr = 1');
    plot(Tb,Yb(:,4), '-','DisplayName','(RK4) Pr = 10');
    plot(Tc,Yc(:,4),'--','DisplayName','(ode45) Pr = 1');
    plot(Td,Yd(:,4),'--','DisplayName','(ode45) Pr = 10');
    xlabel('\eta');
    ylabel('\theta');
    hleg = legend('show');
    set(hleg,'Location','northeastoutside');
    
end

function dy = convection(~, y, Pr)

    %%%%%%
    % Function relating solution variables to their derivatives in the equation for free
    %  convection along a vertical plate.
    %%%

    %  y = F , F' , F'' , theta , theta'
    % dy = F', F'', F''', theta', theta''
    
    dy = zeros(5,1);
    dy(1) = y(2);
    dy(2) = y(3);
    dy(3) = -y(4) + 2 * y(2).^2 - 3 * y(1) .* y(3);
    dy(4) = y(5);
    dy(5) = -3 * Pr * y(1) .* y(5);
end

function val = theta10(Fpp0, Pr, method)

    %%%%%%
    % Returns the value of theta(10) reached by RK4 integration using the given values of
    %  F''(0) and Prandtl number.
    %%%
    
    if Pr == 1
        thp0 = -0.5671;
    elseif Pr == 10
        thp0 = -1.17;
    else
        error('Prandl number must be 1 or 10');
    end
    
    initials = [0, 0, Fpp0, 1, thp0];
    odefun = @(t, y) convection(t, y, Pr);
    if method == 1
        [~,Y] = RK4(odefun, [0,10], 500, initials);
    elseif method == 2
        [~,Y] = ode45(odefun, [0,10], initials);
    else
        error('Method must be 1 (RK4) or 2 (ode45).');
    end
    val = Y(end,2);
end