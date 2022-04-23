function [R, x] = GPT(A, b)
    n = length(A);
    P = [1 : n];
    A = [A, b];
    for p = 1 : n-1
        [_, row] = max(max(A(p : n, p : n), [], 2));
        [_, col] = max(max(A(p : n, p : n), [], 1));
        
        swap = A(p, :);
        A(p, :) = A(p + row - 1, :);
        A(p + row - 1, :) = swap;
        
        swap = A(:, p);
        A(:, p) = A(:, p + col - 1);
        A(:, p + col - 1) = swap;
        
        swap = P(p);
        P(p) = P(p + col - 1);
        P(p + col - 1) = swap;
        
        for i = p+1 : n
            A(i, :) = A(i, :) - A(i, p) / A(p, p) * A(p, :);
        endfor
    endfor
    
    R = A(:, 1 : end-1);
    b = A(:, end);
    x_star = SST(R, b);
    x = zeros(n, 1);
    for p = 1:n
        x(P(p)) = x_star(p);
    endfor
endfunction