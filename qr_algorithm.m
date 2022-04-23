function [X, l, steps] = qr_algorithm(A, max_steps)
    prev_A = A;
    steps = 0;
    Q_final = eye(length(A));
    while 1
        [Q, R] = qr(prev_A);
        A = R * Q;
        Q_final *= Q;
        if (++steps == max_steps)
            l = [];
            for i = 1:length(A)
                if i > 1 && abs(A(i, i-1)) > 1e-4
                    continue
                elseif i < length(A) && abs(A(i + 1, i)) > 1e-4
                    l = [l, eig(A(i:i+1, i:i+1))];
                else
                    l = [l, A(i, i)];
                endif
            endfor
            X = Q_final;
            return;
        endif
        prev_A = A;
    endwhile
endfunction