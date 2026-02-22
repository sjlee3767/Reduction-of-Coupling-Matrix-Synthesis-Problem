function ret = comp_boxes(M1, M2, i, k, eps)
    % i: first box index
    % k: number of boxes

    i_init = i;

    n = size(M1, 1);

    sen1 = 1;
    sen2 = 1;
    for cnt = 1:k
        B1 = M1(i:i+1, i+2:i+3);
        B2 = M2(i:i+1, i+2:i+3);

        D = abs(B1)-abs(B2);
        if max(abs(D), [], 'all') > eps % if two boxes differ, return 0;
            sen1 = 0;
            break;
        end

        i = i+2;
    end

    i = i_init;
    % check for index-reversed symmetry
    for cnt = 1:k
        B1 = M1(i:i+1, i+2:i+3);
        B2 = M2(n-i:n-i+1, n-i-2:n-i-1)';
        tmp = B2(1,1);
        B2(1, 1) = B2(2, 2);
        B2(2, 2) = tmp; % get the box of the opposite direction


        D = abs(B1)-abs(B2);
        if max(abs(D), [], 'all') > eps % if two boxes differ, return 0;
            sen2 = 0;
            break;
        end

        i = i+2;
    end

    if sen1 == 0 && sen2 == 0
        ret = 0;
    else
        ret = 1;
    end
end