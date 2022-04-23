function [r, err, steps] = bisection(f, a, b, tol, max_steps)
    steps = 0;
    while 1
        r = (a + b) / 2;
        if feval(f, a) * feval(f, r) < 0
            b = r;
        else
            a = r;
        endif
        
        err = b - a;
        if (++steps == max_steps) || (err < tol)
            return;
        endif
    endwhile
endfunction