function [R, x] = GPPS(A, b)
    n = length(A);
    A = [A, b];
    for p = 1 : n-1
        s = max(abs(A(p : n, p : n))')';
        [_, ind] = max(abs(A(p : n, p)) ./ s);
        swap = A(p, :);
        A(p, :) = A(p + ind - 1, :);
        A(p + ind - 1, :) = swap;
        
        for i = p+1 : n
            A(i, :) = A(i, :) - A(i, p) / A(p, p) * A(p, :);
        endfor
    endfor
    
    R = A(:, 1 : end-1);
    b = A(:, end);
    x = SST(R, b);
endfunction