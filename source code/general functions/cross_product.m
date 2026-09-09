function ret = cross_product(V)
    n = size(V, 1);

    V = V';

    v = zeros(n, 1);

    for i = 1:n
        if i > 1 && i < n
            subV = [V(:,1:i-1) V(:,i+1:end)];
        elseif i == 1
            subV = V(:, 2:end);
        else
            subV = V(:, 1:n-1);
        end

        v(i) = (-1)^(i-1) * det(subV);
    end

    ret = v;

end