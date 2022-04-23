function [R_] = hessenberg(A)
    [M, N] = size(A);
    R_ = A;
    limit = min(M - 1, N);
    
    for p = 1 : limit - 1
        v = zeros(M, 1);
        v(p + 1 : M, 1) = R_(p + 1 : M, p);
        v(p + 1) += -sign(R_(p + 1, p)) * norm(v);
        
        H = eye(M) - 2 * v * v' / (v' * v);
        
        % Produsul H * R_
        R_ = H * R_ * H';
    end
end