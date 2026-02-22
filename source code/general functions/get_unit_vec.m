function ret = get_unit_vec(th)
    n = length(th);
    ret = ones(n+1, 1);
    for i = 1:n
        ret(i) = ret(i) * cos(th(i));
        for j = i+1:n+1
            ret(j) = ret(j) * sin(th(i));
        end
    end
end