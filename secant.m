function [r, err, steps] = secant(f, x0, x1, tol, max_steps)
    steps = 0;
    prev_r = x1;
    prev_prev_r = x0;
    while 1
        r = prev_r - (prev_prev_r - prev_r) / (f(prev_prev_r) - f(prev_r)) * f(prev_r);
        
        err = abs(r - prev_r);
        abs(sqrt(3) - r)
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_prev_r = prev_r;
        prev_r = r;
    endwhile
endfunction