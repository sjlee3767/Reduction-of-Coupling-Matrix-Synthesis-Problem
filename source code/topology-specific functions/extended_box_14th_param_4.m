function [optimizing_param, M, R] = extended_box_14th_param_4(M_init, th)
    % 5 tzs

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

    M = M_init;

    N = length(M);
    
    I = eye(N);
    
    P = I;
    
    r0 = I(:, 1);
    r1 = I(:, 2);
    r8 = I(:, end-1);
    r9 = I(:, end);

    r7 = zeros(10, 1);
    r7(4:8) = get_unit_vec_spherical(th);

    V = [r0 r1 r7 r8 r9];
    r6 = normalize_C(get_next_orthdcp([V], M*r8));
    V = [V, r6];
    r5 = normalize_C(get_next_orthdcp([V ], M*r7));
    V = [V, r5];
    r4 = normalize_C(get_next_orthdcp([V ], M*r6));
    V = [V, r4];
    r3 = normalize_C(get_next_orthdcp([V ], M*r5));
    V = [V, r3];
    r2 = normalize_C(get_next_orthdcp([V ], M*r4));
    
    R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9];
    M = R'*M*R;

    M = M_trim_zeros(M);

    % tmp = (M.*(ones(N)-T_target)).^2;
    % optimizing_param = sqrt(sum(tmp(:)));


    optimizing_param = [];
    for i = 1:N
        for j = i:N
            if T(i, j) == 0
                optimizing_param = [optimizing_param, M(i, j)^2];
            end
        end
    end
    
    
end
