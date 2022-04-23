function [x, err, steps] = matrix_SOR(A, b, x0, tol, max_steps, w)
    N = diag(diag(A)) + w * tril(A, -1);
    P = (1 - w) * diag(diag(A)) - w * triu(A, 1);
    inv_N = inv(N);
    G = inv_N * P;
    c = w * inv_N * b;
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