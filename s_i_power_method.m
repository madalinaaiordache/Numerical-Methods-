function [lambda, x, err, steps] = s_i_power_method(A, x0, tol, max_steps, alfa)
    prev_x = x0;
    steps = 0;
    
    while 1
        x = (A - alfa * eye(length(A))) \ prev_x;
        x = x / norm(x);
        lambda = x' * A * x;
    
        err = min(norm(x - prev_x), norm(x + prev_x));
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_x = x;
    endwhile
endfunction