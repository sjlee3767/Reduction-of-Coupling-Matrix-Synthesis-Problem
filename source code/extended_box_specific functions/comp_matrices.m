function ret = comp_matrices(M1, M2, eps)
    n = size(M1, 1);

    sen1 = 1;
    sen2 = 1;


    D = abs(M1)-abs(M2);

    if max(abs(D), [], 'all') > eps % if two matrices differ, return 0;
        sen1 = 0;
    end

    
    M3 = M2; % M3: index-reversed version of M2.
    for i = 1:n
        for j = 1:n
            M3(i, j) = M2(n-i+1, n-j+1);
        end
    end

    D = abs(M1)-abs(M3);

    if max(abs(D), [], 'all') > eps % if two matrices differ, return 0;
        sen2 = 0;
    end

    if sen1 == 0 && sen2 == 0
        ret = 0;
    else
        ret = 1;
    end
end