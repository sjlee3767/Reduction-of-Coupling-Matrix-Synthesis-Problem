function T = proj(A)
    m = size(A, 2);
    
    if m < 1
        T = eye(size(A, 1));
        return;
    end
    T = A(:, 1)*A(:, 1)';
    for i = 2:m
        tmp = A(:, i);
        T = T + tmp*tmp';
    end

end