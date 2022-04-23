function [x] = SIT(A, b)
    n = length(b);
    x = [b(1) / A(1, 1)];
    for i = 2:n
        x = [x; (b(i) - A(i, [1:i-1]) * x) / A(i,i)];
    endfor
endfunction