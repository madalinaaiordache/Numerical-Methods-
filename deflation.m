function [X, l, total_steps] = deflation(A, tol, max_steps)
    X = [];
    l = [];
    A0 = A;
    total_steps = 0;
    n = length(A);
    for i = 1:n
        [lambda, x_B, err, steps] = power_method(A, rand(n - i + 1, 1), tol, max_steps);
        total_steps += steps;
        [lambda, x, err, steps] = s_i_power_method(A0, rand(n, 1), tol, max_steps, lambda + rand() * 1e-2);
        l = [l, lambda];
        X = [X, x];
        if i < n
            v = x_B; v(1) += -sign(x_B(1)) * norm(x_B);
            H = eye(n - i + 1) - 2 * (v * v') / (v' * v);
            A = H * A * H';
            A = A(2:end, 2:end);
        endif
    endfor
endfunction