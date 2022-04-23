function [X, l, total_steps] = orthogonal_power_method(A, tol, max_steps)
    X = [];
    l = [];
    total_steps = 0;
    n = length(A);

    for i = 1:n
        prev_x = rand(n, 1);
        steps = 0;
        while 1
            % Orthogonalize span{X(:, 1), X(:, 2), ...} out of prev_x
            for j = 1:size(X, 2)
                prev_x -= X(:, j)' * prev_x * X(:, j);
            endfor
            
            % Power method
            x = A * prev_x;
            x = x / norm(x);
            lambda = x' * A * x;

            err = norm(x - prev_x);
            if (++steps == max_steps) || (err < tol)
                l = [l, lambda];
                X = [X, x];
                total_steps += steps;
                break;
            endif

            prev_x = x;
        endwhile
    endfor
endfunction