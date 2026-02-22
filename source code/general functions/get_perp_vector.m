% gets a set of column vectors A and a vector u.
% 
function [v, A] = get_perp_vector(A, u)
    n = size(A, 1);
    if n == 0
        n = length(u);
        T = eye(n);
    else
        I = eye(n);
        T = I - proj(A);
    end
    
    
    v = T*u;
    v = normalize_C(v);
    A = [A, v];
end