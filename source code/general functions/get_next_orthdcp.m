function ret = get_next_orthdcp(V, v, opt) % returns a vector orthogonal to all column vectors of V so that they span the vector v. (column vectors of V are orthogonal to each other)
    [n, m] = size(V);
    if n == 0 && m == 0
        ret = v;
        return;
    end
    if nargin <= 2
        opt = 0;
    end
    if opt == 0
        if m >= n
            disp("warning! number of vectors exceed the dimension!");
            tmp = rank(V);
            if tmp >= n
                disp("the set of vectors to be orthogonal span the whole space!");
                % keyboard;
                return;
            end
        end
    end

    if n ~= length(v)
        disp("dimensions of vectors do not match!");
        return;
    end

    c = zeros(m, 1);
    for i = 1:m
        if V(:, i)'*V(:, i) < 1e-10
            continue;
        end
        c(i) = (V(:, i)' * v) / (V(:, i)'*V(:, i));
    end
    ret = v - V*c;
end

