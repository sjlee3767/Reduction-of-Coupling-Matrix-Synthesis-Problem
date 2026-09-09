function [optimizing_param, M, R, depth] = appendix_ex_1(M_init, th)
% gets 2 variables

    M = M_init;

    N = length(M);

    I = eye(N);

    P = I;

    % r0 = I(:, 1);
    % r11 = I(:, end);

    r0 = I(:, 1);
    r9 = I(:, end);
    V = [r0 r9];
    r1 = normalize_C(get_next_orthdcp(V, M*r0));
    r8 = normalize_C(get_next_orthdcp(V, M*r9));
    V = [V r1 r8];
    A = M(3:end-2,3:end-2);
    [v, lam] = eig_C(A);

    r4 = zeros(N, 1);
    global branch;
    r4(3:end-2) = v(:,branch);
    V = [V r4];

    r2 = normalize_C(get_next_orthdcp(V, M*r1));
    V = [V r2];
    r3 = normalize_C(get_next_orthdcp(V, M*r2));
    V = [V r3];

    P1 = align_vector(r2, [0 0 2 1 1 1 1 1 0 0]);
    P = P1*P;
    P2 = align_vector(P*r3, [0 0 0 2 1 1 1 1 0 0]);
    P = P2*P;
    P3 = align_vector(P*r4, [0 0 0 0 2 1 1 1 0 0]);
    P = P3*P;
    Pr5 = zeros(N, 1);
    global initial_point;
    global tangent;
    tmp = initial_point' + tangent*th;
    Pr5(6:8) = normalize_C(tmp);
    r5 = P'*Pr5;
    V = [V r5];

    r6 = normalize_C(get_next_orthdcp(V, M*r5));
    V = [V r6];

    r7 = cross_product(V);

    R = [ r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 ];

    M = R'*M*R;

    M = M_trim_zeros(M);
    % return;
    % M = M_annihilate2(M, 4, 7, 4, 8);
    % disp_topology(M);
    % M = M_annihilate2(M, 4, 6, 4, 7);
    % disp_topology(M);
    % M = M_annihilate2(M, 3, 6, 3, 7);
    % disp_topology(M);
    % M = M_annihilate2(M, 3, 4, 3, 6);
    % disp_topology(M);

    actions = [4, 7, 4, 8, 1
                   4, 6, 4, 7, 1
                   3, 6, 3, 7, 1
                   3, 4, 3, 6, 1];

    tmp = bst2(M, actions, 0, th, 1);
    M = tmp{1};
    optimizing_param = tmp{2};
    depth = tmp{3};

    M = M_trim_zeros(M);

    % M = M_init;
    % 
    % N = length(M);
    % 
    % I = eye(N);
    % 
    % P = I;
    % 
    % r0 = I(:, 1);
    % r11 = I(:, end);
    % 
    % r1 = I(:, 2);
    % r10 = I(:, end-1);
    % V = [r0 r1 r10 r11];
    % r2 = normalize_C(get_next_orthdcp(V, M*r1));
    % r9 = normalize_C(get_next_orthdcp(V, M*r10));
    % V = [V r2 r9];
    % A = M(4:end-3,4:end-3);
    % [v, lam] = eig_C(A);
    % 
    % r5 = zeros(N, 1);
    % global branch;
    % r5(4:end-3) = v(:,branch);
    % V = [V r5];
    % 
    % r3 = normalize_C(get_next_orthdcp(V, M*r2));
    % V = [V r3];
    % r4 = normalize_C(get_next_orthdcp(V, M*r3));
    % V = [V r4];
    % 
    % P1 = align_vector(r3, [0 0 0 2 1 1 1 1 1 0 0 0]);
    % P = P1*P;
    % P2 = align_vector(P*r4, [0 0 0 0 2 1 1 1 1 0 0 0]);
    % P = P2*P;
    % P3 = align_vector(P*r5, [0 0 0 0 0 2 1 1 1 0 0 0]);
    % P = P3*P;
    % 
    % Pr6 = zeros(N, 1);
    % global initial_point;
    % global tangent;
    % tmp = initial_point' + tangent*th;
    % Pr6(7:9) = normalize_C(tmp);
    % r6 = P'*Pr6;
    % V = [V r6];
    % 
    % r7 = normalize_C(get_next_orthdcp(V, M*r6));
    % V = [V r7];
    % 
    % r8 = cross_product(V);
    % 
    % R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11];
    % 
    % M = R'*M*R;
    % 
    % M = M_trim_zeros(M);
    % 
    % % M = M_annihilate2(M, 4, 7, 4, 8);
    % % disp_topology(M);
    % % M = M_annihilate2(M, 4, 6, 4, 7);
    % % disp_topology(M);
    % % M = M_annihilate2(M, 3, 6, 3, 7);
    % % disp_topology(M);
    % % M = M_annihilate2(M, 3, 4, 3, 6);
    % % disp_topology(M);
    % 
    % actions = [4, 7, 4, 8, 1
    %                4, 6, 4, 7, 1
    %                3, 6, 3, 7, 1
    %                3, 4, 3, 6, 1];
    % actions(:, 1:end-1) = actions(:, 1:end-1)+1; % change to 1-based indexing
    % 
    % tmp = bst2(M, actions, 0, th, 1);
    % M = tmp{1};
    % optimizing_param = tmp{2};
    % depth = tmp{3};
    % 
    % M = M_trim_zeros(M);
    
end
    