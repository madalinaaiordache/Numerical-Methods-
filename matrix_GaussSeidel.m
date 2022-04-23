function [x, err, steps] = matrix_GaussSeidel(A, b, x0, tol, max_steps)
    N = tril(A);
    P = N - A;
    inv_N = inv(N);
    G = inv_N * P;
    c = inv_N * b;
    prev_x = x0;
    
    steps = 0;
    while 1
        x = G * prev_x + c;
        
        err = norm(x - prev_x);
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_x = x;
    endwhile
endfunction