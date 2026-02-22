function ret = get_unit_vec_spherical(th)
    n = length(th)+1;
    v = ones(n, 1);
    for i = 1:n-1
        for j = 1:n-i
            v(j) = v(j)*cos(th(i));
        end
        v(n-i+1) = v(n-i+1)*sin(th(i));
    end

    ret = v;
end