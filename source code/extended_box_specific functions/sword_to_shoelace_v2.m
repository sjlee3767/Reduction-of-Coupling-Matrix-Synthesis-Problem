function ret = sword_to_shoelace_v2(M)
    sympref('FloatingPointOutput',true);

    % sympref('default');

    M_init = M;
    
    N = length(M);
    n = N-2;
    
    I = eye(N);

    l = 7;
    r = n/2+3;
    k = (r-l)-5;

    p = 3;
    while r < n-1
        M = M_annihilate(M, p, p+1, p, p+1);
        M = level_1(M, l+1, r+1, k);
        

        l = l + 2;
        r = r + 1;
        k = k - 1;
        p = p+2;
    end

    % ret = M;
    % return;
    for i = p:2:N-3
        M = M_annihilate(M, i, i+1, i, i+1);
    end

    ret = M;


end

function ret = level_1(M, l, r, k) % k+1: number of level_2 to be performed
    for i = 0:k
        M = level_2(M, r, r-i-1);
    end
    
    for i = 1:r-l % r-l: number of column elements to be annihilated
        if (r)-(l+i-1) > 3
            M = M_annihilate(M, l+1+i, l+2+i, l+i, l+2+i);
        end
        
        
        M = M_annihilate(M, l+i-1, r, l-3+i, r);
        
        if (r)-(l+i-1) > 3
            M = M_annihilate(M, l+1+i, l+2+i, l-1+i, l+2+i);
        end
        
        if i < k+1 % do level_2_inv only when level_2 was actually performed
            M = level_2_inv(M, r, (r-k-1)+i-1);
        end
    end

    ret = M;

end

function ret = level_2(M, r, k) % k: index of row to be annihilated
    M = M_annihilate(M, k-1, k, k, r);
    for i = 1:r-k-2
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, k, r-i);
        M = level_3_inv(M, r, i);
    end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, r-i-2, r-i);
        M = level_3_inv(M, r, i);
     end

     ret = M;
end

function ret = level_2_inv(M, r, k) % k: index of row to be revived
     for i = 1:r-k-2
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, k-2, r-i);
        M = level_3_inv(M, r, i); 
     end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, r-i-2, r-i);
        M = level_3_inv(M, r, i);
     end

     M = M_annihilate(M, k-1, k, k-3, k);

     ret = M;
end

function ret = level_3(M, r, k)
    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = 3:2:k
            M = M_annihilate(M, r-i+1, r-i+2, r-i+1, r-i+2);
        end
    else
        for i = 4:2:k
            M = M_annihilate(M, r-i+1, r-i+2, r-i+1, r-i+2);
        end
    end
    ret = M;
end

function ret = level_3_inv(M, r, k)
    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = k:-2:3
            M = M_annihilate(M, r-i+1, r-i+2, r-i-1, r-i+2);
        end
    else
        for i = k:-2:4
            M = M_annihilate(M, r-i+1, r-i+2, r-i-1, r-i+2);
        end
    end
    ret = M;
end