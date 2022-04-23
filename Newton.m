function [R, err, steps] = secant(F_handlers, X0, tol, max_steps)
    steps = 0;
    prev_R = X0;
    DF = zeros(length(F_handlers), length(X0));
    F = zeros(length(F_handlers), 1);
    h = 1e-5;
    
    while 1
        for i = 1:size(DF, 1)
            for j = 1:size(DF, 2)
                R_plus = prev_R; R_plus(j) += h;
                R_minus = prev_R; R_minus(j) -= h;
                DF(i, j) = 1 / (2*h) * (F_handlers{i}(R_plus) - F_handlers{i}(R_minus));
            endfor
        endfor
        
        for i = 1:length(F)
            F(i) = F_handlers{i}(prev_R);
        endfor
        
        R = prev_R - DF \ F;
        err = norm(R - prev_R);
        if (++steps == max_steps) || (err < tol)
            return;
        endif
        
        prev_R = R;
    endwhile
endfunction