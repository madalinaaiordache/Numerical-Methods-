function [L, U] = Crout(A)
    n = length(A);
    L = zeros(n);
    U = eye(n);
    for p = 1:n
        for i = p:n
            L(i, p) = A(i, p) - L(i, [1:p-1]) * U([1:p-1], p);
        endfor
        for j = p+1:n
            U(p, j) = 1 / L(p, p) * (A(p, j) - L(p, [1:p-1]) * U([1:p-1], j));
        endfor
    endfor
endfunction