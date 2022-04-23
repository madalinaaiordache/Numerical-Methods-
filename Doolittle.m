function [L, U] = Doolittle(A)
    n = length(A);
    L = eye(n);
    U = zeros(n);
    for p = 1:n
        for j = p:n
            U(p, j) = A(p, j) - L(p, [1:p-1]) * U([1:p-1], j);
        endfor
        for i = p+1:n
            L(i, p) = 1 / U(p, p) * (A(i, p) - L(i, [1:p-1]) * U([1:p-1], p));
        endfor
    endfor
endfunction