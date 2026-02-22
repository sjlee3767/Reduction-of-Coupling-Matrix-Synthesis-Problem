function [optimizing_param, M, R, depth] = extended_box_8th_C(M_init, th)
    % 3 tzs
    M = M_init;

    N = length(M);
    
    I = eye(N);
    
    P = I;
    
    r0 = I(:, 1);
    r1 = I(:, 2);
    r2 = I(:, 3);
    r7 = I(:, end-2);
    r8 = I(:, end-1);
    r9 = I(:, end);
    
    P1 = align_vector(M*r2, [0 0 0 2 1 1 0 0 0 0]);
    P = P1*P;
    
    Pr6 = zeros(10, 1);
    Pr6(5:7) = get_unit_vec_spherical(th);
    
    r6 = P'*Pr6;
    get_next_orthdcp([r6, r7, r8], M*r7);
    r5 = normalize_C(get_next_orthdcp([r6, r7, r8], M*r7));
    r4 = normalize_C(get_next_orthdcp([r5, r6, r7], M*r6));
    r3 = normalize_C(get_next_orthdcp([r2, r4, r5, r6, r7], M*r5));
    
    R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9];
    M = R'*M*R;

    % M = M_annihilate(M, 3, 4, 3, 4, 0);
    % M = M_annihilate(M, 5, 6, 5, 6, 0);
    % M = M_annihilate(M, 7, 8, 7, 8, 0);

    actions = [3, 4, 3, 4, 1
        5, 6, 5, 6, 1
        7, 8, 7, 8, 1];

    tmp = bst(M, actions, 0, 2, th, 1);
    M = tmp{1};
    optimizing_param = tmp{2};
    depth = tmp{3};
    
    M = M_trim_zeros(M);

    % optimizing_param = [get_box_absmin(M, 3), get_box_absmin(M, 5)];
end
