function [L] = Cholesky(A)
    n = length(A);
    L = NaN;
    
    if !isequal(A, A')
        return;
    endif
    
    for i = 1:n
        if det(A(1:i, 1:i)) <= 0
            return;
        endif
    endfor
    
    L = zeros(n);
    for i = 1:n
        for j = 1:i-1
            L(i, j) = (A(i, j) - L(i, [1:j-1]) * L(j, [1:j-1])') / L(j, j);
        endfor
        L(i, i) = sqrt(A(i, i) - L(i, [1:i-1]) * L(i, [1:i-1])');
    endfor
endfunction