function [x, err, steps] = Jacobi(A, b, x0, tol, max_steps)
    prev_x = x0;
    x = prev_x;
    n = length(x);
    
    steps = 0;
    while 1
        for i = 1 : n
            x(i) = (b(i) - A(i, :) * prev_x + A(i, i) * prev_x(i)) / A(i, i);
        endfor
        
        err = norm(x - prev_x);
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_x = x;
    endwhile
endfunction