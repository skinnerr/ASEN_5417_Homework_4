function [ l, u ] = LU_Decompose( a, b, c )
    
    %%%%%%
    % Decomposes a tri-diagonal matrix equation into L and U matrices.
    %     a -- diagonal
    %     b -- sub-diagonal
    %     c -- super-diagonal
    %     l -- diagonal of L matrix
    %     u -- super-diagonal of U matrix
    %
    % Ryan Skinner, October 2015
    %%%
    
    N = length(a);
    
    l = nan(N,  1);
    u = nan(N-1,1);
    
    l(1) = a(1);
    for i = 2:N
        u(i-1) = c(i-1) / l(i-1);
        l(i) = a(i) - b(i-1) * u(i-1);
    end

end