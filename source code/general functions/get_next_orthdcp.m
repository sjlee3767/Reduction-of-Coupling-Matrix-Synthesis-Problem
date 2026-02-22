function ret = get_next_orthdcp(V, v) % returns a vector orthogonal to all column vectors of V so that they span the vector v. (column vectors of V are orthogonal to each other)
    [n, m] = size(V);
    if n == 0 && m == 0
        ret = v;
        return;
    end
    if m >= n
        disp("number of vectors exceed the demension!");
        return;
    end

    if n ~= length(v)
        disp("dimensions of vectors do not match!");
        return;
    end

    c = zeros(m, 1);
    for i = 1:m
        c(i) = (V(:, i)' * v) / (V(:, i)'*V(:, i));
    end
    ret = v - V*c;
end