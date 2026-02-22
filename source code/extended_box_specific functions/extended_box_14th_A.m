function [optimizing_param, M, R] = extended_box_14th_A(M_init, th)
    M = M_init;

    N = length(M);
    n = N-2;
    
    T = [0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         1     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0
         0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0
         0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0
         0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0
         0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0
         0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0
         0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0
         0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0
         0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0
         0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0
         0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0
         0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0
         0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     0
         0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0];

    I = eye(N);
    
    P = I;
    
    l = 1;
    for i = 2:N-2
        cur_th = th(l:l+n-i);
        Pr_i = zeros(N, 1);
        Pr_i(i:N-1) = get_unit_vec_spherical(cur_th);

        idx_vec = zeros(1, n);
        idx_vec(i-1) = 2;
        for j = i:n
            idx_vec(j) = 1;
        end
        P_tmp = I;
        P_tmp(2:N-1, 2:N-1) = align_vector(Pr_i(2:end-1), idx_vec);

        P = P_tmp*P;

        

        l = l + n-i+1;
    end

    R = P';

    M = R'*M*R;
    

    optimizing_param = [];
    for i = 1:N
        for j = i:N
            if T(i, j) == 0
                optimizing_param = [optimizing_param, M(i, j)^2];
            end
        end
    end
end
