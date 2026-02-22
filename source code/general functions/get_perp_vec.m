function ret = get_perp_vec(V) % returns a vector perpenticular to all n-1 column vectors of V
    n = size(V, 1);
    if size(V, 2) ~= n-1
        disp("dimension does not match number of vectors");
        return;
    end
    ret = zeros(n, 1);
    M = zeros(n-1, n-1);
    V = V';
    for i = 1:n
        M(:, 1:i-1) = V(:, 1:i-1);
        M(:, i:end) = V(:, i+1:end);
        ret(i) = (-1)^(i+1) * det(M);
    end


end