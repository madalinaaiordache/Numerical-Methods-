function [R, x] = GPP(A, b)
    n = length(A);
    A = [A, b];
    for p = 1 : n-1
        [_, ind] = max(abs(A(p, p : n)));
        swap = A(p, :);
        A(p, :) = A(ind, :);
        A(ind, :) = swap;
        
        for i = p+1 : n
            A(i, :) = A(i, :) - A(i, p) / A(p, p) * A(p, :);
        endfor
    endfor
    
    R = A(:, 1 : end-1);
    b = A(:, end);
    x = SST(R, b);
endfunction