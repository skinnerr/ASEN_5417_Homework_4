function [ x ] = Thomas( a, b, c, rhs )
    
    %%%%%%
    % Solves a tri-diagonal matrix system using the Thomas algorithm.
    %     a -- diagonal
    %     b -- sub-diagonal
    %     c -- super-diagonal
    %   rhs -- right-hand side vector
    %   sol -- solution vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    % Eliminate the sub-diagonal.
    for i = 1:length(c)
        r = b(i) / a(i);
          a(i+1) =   a(i+1) - r * c(i);
        rhs(i+1) = rhs(i+1) - r * rhs(i);
    end
    
    % Back-substitute and calculate the solution vector.
    x = zeros(length(a),1);
    x(end) = rhs(end) / a(end);
    for i = length(a)-1:-1:1
        x(i) = (rhs(i) - c(i) * x(i+1)) / a(i);
    end

end