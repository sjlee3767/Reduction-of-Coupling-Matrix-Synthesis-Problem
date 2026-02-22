function [optimizing_param, M, R, depth] = extended_box_14th_C(M_init, th)
    M = M_init;

    N = length(M);
    
    I = eye(N);
    
    P = I;

    r0 = I(:, 1);
    r1 = I(:, 2);
    r2 = I(:, 3);
    r11 = I(:, end-4);
    r12 = I(:, end-3);
    r13 = I(:, end-2);
    r14 = I(:, end-1);
    r15 = I(:, end);

    r10 = normalize_C(get_next_orthdcp([r11 r12], M*r11));
    
    P1 = align_vector(r10, [0 0 0 0 1 1 0 0 0 0 2 0 0 0 0 0]);
    P = P1*P;
    P2 = align_vector(P*M*r10, [0 0 0 1 1 1 1 0 0 2 0 0 0 0 0 0]);
    P = P2*P;

    Pr3 = zeros(N, 1);
    Pr3(4:9) = get_unit_vec_spherical(th);

    r3 = P'*Pr3;

    v1 = normalize_C((I-r1*r1'-r2*r2'-r3*r3')*M*r2);
    v2 = normalize_C((I-r2*r2'-r3*r3')*M*r3);

    v3 = normalize_C((I-v2*v2')*v1); % {v2, v3} is an orthogonal basis of the linear space of {v1, v2}
    Mr10 = M*r10;

    r5 = normalize_C((Mr10'*v3)*v3 + (Mr10'*v2)*v2);
    r4 = normalize_C(get_next_orthdcp([r1 r2 r3 r5], M*r2));

    r6 = normalize_C(get_next_orthdcp([r2 r3 r4 r5], M*r4));
    r7 = normalize_C(get_next_orthdcp([r2 r3 r4 r5 r6 r10], M*r5));
    r8 = normalize_C(get_next_orthdcp([r4 r5 r6 r7 r10], M*r6));
    r9 = normalize_C(get_next_orthdcp([r5 r6 r7 r8 r10], M*r7));
    
    
    
    R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15];
    
    M = R'*M*R;

    t1 = toc;
    
    % M = M_annihilate(M, 8, 11, 6, 11);
    % M = M_annihilate(M, 9, 11, 7, 11);
    % M = M_annihilate(M, 10, 11, 8, 11);
    % M = M_annihilate(M, 10, 12, 8, 12);
    % M = M_annihilate(M, 11, 12, 9, 12);
    % M = M_annihilate(M, 12, 13, 10, 13);
    % 
    % M = M_annihilate(M, 5, 6, 5, 6, 0);
    % M = M_annihilate(M, 7, 8, 7, 8, 0);
    % M = M_annihilate(M, 9, 10, 9, 10, 0);
    % M = M_annihilate(M, 11, 12, 11, 12, 0);
    % M = M_annihilate(M, 13, 14, 13, 14, 0);
    % M = M_annihilate(M, 3, 4, 3, 4, 0);

    actions = [8, 11, 6, 11, 1
        9, 11, 7, 11, 1
        10, 11, 8, 11, 1
        10, 12, 8, 12, 1
        11, 12, 9, 12, 1
        12, 13, 10, 13, 1
        5, 6, 5, 6, 1
        7, 8, 7, 8, 1
        9, 10, 9, 10, 1
        11, 12, 11, 12, 1
        13, 14, 13, 14, 1
        3, 4, 3, 4, 1];

    tmp = bst(M, actions, 0, 5, th, 1);
    M = tmp{1};
    optimizing_param = tmp{2};
    depth = tmp{3};
    
    
    % M = M_annihilate(M, 8, 9, 8, 10)
    % M = M_annihilate(M, 6, 7, 6, 9)
    % M = M_annihilate(M, 4, 5, 4, 7)
    % M = M_annihilate(M, 2, 3, 1, 3)

    
    M = M_trim_zeros(M);

    % optimizing_param = zeros(1, 5);
    % for i = 1:5
    %     optimizing_param(i) = get_box_absmin(M, i*2+1);
    % end

    % t2 = toc;
    % keyboard;
end






