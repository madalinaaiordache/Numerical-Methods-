function [X, l, steps] = orthogonal_iteration(A, tol, max_steps)
    prev_X = rand(length(A));
    steps = 0;
    while 1
        [Q, R] = qr(prev_X);
        X = A * Q;
        if (norm(X - prev_X) < tol) || (++steps == max_steps)
            [Q, R] = qr(X);
            X = Q;
            l = diag(Q' * A * Q)';
            return;
        endif
        prev_X = X;
    endwhile
endfunction