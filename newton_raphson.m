function [r, err, steps] = newton_raphson(f, df, x0, tol, max_steps)
    steps = 0;
    prev_r = x0;
    while 1
        r = prev_r - f(prev_r) / df(prev_r)
        
        err = abs(r - prev_r);
        abs(sqrt(3) - r)
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_r = r;
    endwhile
endfunction