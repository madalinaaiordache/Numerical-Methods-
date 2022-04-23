function [res] = custom_norm(x, p)
    if p == inf
        res = max(abs(x));
    else
        res = sum(abs(x) .^ p)^(1/p);
    endif
endfunction