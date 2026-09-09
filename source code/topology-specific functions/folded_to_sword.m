function ret = folded_to_sword(M, th)
    M_init = M;

    N = length(M);
    n = N-2;
    
    % th = rand(1, n/2-2) * 2*pi;
    
    if length(th) ~= n/2-2 % = 2+(n-8)/2
        fprintf("Number of free variables (ths) does not match the size of the matrix M!\n");
        ret = 0;
        return;
    end
    
    I = eye(N);
    R = zeros(N);
    A = [];
    
    R(:, 1) = I(:, 1);
    R(:, N) = I(:, N);
    
    A = [R(:, 1), R(:, N)];
    
    for k = 0:1
        [R(:, k+2), A] = get_perp_vector(A, M*R(:, k+1));
    end
    
    for k = N-1:-1:N/2+3
        [R(:, k), A] = get_perp_vector(A, M*R(:, k+1));
    end
    
    P = I;
    % idx_vec = [zeros(1, N/2-2), ones(1, N/2+2)];
    idx_vec = [ones(1, N-1), 2];
    % keyboard;
    B = [A, M*R(:,N/2+3)];
    for i = 1:size(B, 2)
        % keyboard;
        P_tmp = align_vector(P*B(:, i), idx_vec);
        P = P_tmp*P;
        idx_vec(N-i+1) = 0;
        idx_vec(N-i) = 2;
    
        % idx_vec(N/2-2+i-1) = 0;
        % idx_vec(N/2-2+i) = 2;
        % keyboard;
        % P_tmp = align_vector(P*B(:, i), idx_vec);
        % P = P_tmp*P;
    end
    
    Pr3 = zeros(N, 1);
    Pr3(1:N/2-2) = get_unit_vec_spherical(th);
    R(:,4) = P'*Pr3;
    A = [A, R(:,4)];
    
    % return;
    
    v1 = get_perp_vector([R(:, 2), R(:, 3), R(:, 4)], M*R(:, 3));
    v2 = get_perp_vector([R(:, 3), R(:, 4)], M*R(:, 4));
    v3 = get_perp_vector([v2], v1);
    R(:, 6) = normalize_C(proj([v2,v3])*M*R(:,N/2+3));
    A = [A, R(:, 6)];
    [R(:, 5), A] = get_perp_vector(A, M*R(:, 3));
    
    for k=4:N/2-1
        [R(:,k+3), A] = get_perp_vector(A, M*R(:, k+1));
    end
    
    % array2table(M_trim_zeros(R'*R))
    
    
    ret = R'*M*R;
    
    % array2table(M_trim_zeros(ret))

end