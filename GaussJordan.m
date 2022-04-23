function [invA] = GaussJordan(A)
    n = length(A);
    A = [A, eye(n)];
    
    for p = 1 : n
        A(p, :) = A(p, :) / A(p, p);
        
        for i = 1 : n
            if i != p
                A(i, :) = A(i, :) - A(i, p) / A(p, p) * A(p, :);
            endif
        endfor
    endfor
    
    invA = A(:, n+1 : 2*n);
endfunction