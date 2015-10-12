function [ F ] = Euler( g, F0, delta )
    
    %%%%%%
    % Integrates the data set g using the Euler method, given the initial value F0 and the
    % uniform spacing delta between grid points.
    %
    % Ryan Skinner, October 2015
    %%%
    
    F = nan(length(g),1);
    F(1) = F0;
    for i = 1:(length(g)-1)
        F(i+1) = F(i) + delta * g(i);
    end
    
end