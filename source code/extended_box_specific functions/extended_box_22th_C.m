function [optimizing_param, M, R, depth] = extended_box_22th_C(M_init, th) % 9 variables
    % tic;
    M = M_init;
    
    N = length(M);
    
    I = eye(N);
    
    
    P = I;

    r0 = I(:, 1);
    r1 = I(:, 2);
    r2 = I(:, 3);
    r17 = I(:, end-6);
    r18 = I(:, end-5);
    r19 = I(:, end-4);
    r20 = I(:, end-3);
    r21 = I(:, end-2);
    r22 = I(:, end-1);
    r23 = I(:, end);
    
    
    r16 = normalize_C(get_next_orthdcp([r17 r18],  M*r17));
    r15 = normalize_C(get_next_orthdcp([r16 r17], M*r16));
    r14 = normalize_C(get_next_orthdcp([r15 r16], M*r15));
    Mr14 = M*r14;
    
    P1 = align_vector(r16,      [0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0]);
    P = P1*P;
    P2 = align_vector(P*r15,    [0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0]);
    P = P2*P;
    P3 = align_vector(P*r14,    [0 0 0 0 1 1 1 1 1 1 0 0 0 0 2 0 0 0 0 0 0 0 0 0]);
    P = P3*P;
    P4 = align_vector(P*Mr14,   [0 0 0 1 1 1 1 1 1 1 1 0 0 2 0 0 0 0 0 0 0 0 0 0]);
    P = P4*P;
    
    Pr3 = zeros(N, 1);
    Pr3(4:13) = get_unit_vec_spherical(th);
    
    r3 = P'*Pr3;
    
    v1 = normalize_C(get_next_orthdcp([r1 r2 r3], M*r2));
    v2 = normalize_C(get_next_orthdcp([r2 r3], M*r3));
    
    v3 = normalize_C(get_next_orthdcp([v2], v1)); % {v2, v3} is an orthogonal basis of the linear space of {v1, v2}
    
    r5 = normalize_C((Mr14'*v3)*v3 + (Mr14'*v2)*v2);
    r4 = normalize_C(get_next_orthdcp([r1 r2 r3 r5], M*r2));
    
    
    r6 = normalize_C(get_next_orthdcp([r2 r3 r4 r5], M*r4));
    r7 = normalize_C(get_next_orthdcp([r2 r3 r4 r5 r6 r14], M*r5));
    r8 = normalize_C(get_next_orthdcp([r4 r5 r6 r7 r14], M*r6));
    r9 = normalize_C(get_next_orthdcp([r5 r6 r7 r8 r14], M*r7));
    r10 = normalize_C(get_next_orthdcp([r6 r7 r8 r9 r14], M*r8));
    r11 = normalize_C(get_next_orthdcp([r7 r8 r9 r10 r14], M*r9));
    r12 = normalize_C(get_next_orthdcp([r8 r9 r10 r11 r14], M*r10));
    r13 = normalize_C(get_next_orthdcp([r9 r10 r11 r12 r14], M*r11));
    
    
    
    R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19 r20 r21 r22 r23];
    
    M = R'*M*R;
    % t1 = toc;

    % keyboard;
    % transform to shoelace
    % k = 3;
    % for idx = 0:6
    %     i = idx+N/2+3;
    %     d = idx;
    %     for j = 1:k
    %         M = M_annihilate(M, i-1-j, i-j, i-j, i);
    %     end
    % 
    % 
    %     p = 8+d;
    % 
    %     for j = p:p+6
    %         if k == 3 && j == p
    %             M = M_annihilate(M, i-2, i-1, i-3, i-1);
    %         end
    %         if i-j > 3
    %             M = M_annihilate(M, j+2, j+3, j+1, j+3);
    %         end
    %         M = M_annihilate(M, j, i, j-2, i);
    %         if i-j > 3
    %             M = M_annihilate(M, j+2, j+3, j, j+3);
    %         end
    %         if k == 3 && j == p
    %             M = M_annihilate(M, i-2, i-1, i-5, i-1);
    %         end
    %     end
    %     k = k-1;
    % end
    % 
    % for i = 3:2:N-3
    %     M = M_annihilate(M, i, i+1, i, i+1);
    % end
    % 
    % optimizing_param = zeros(1, 9);
    % for i = 1:9
    %     optimizing_param(i) = get_box_absmin(M, i*2+1)^2;
    % end
    % depth = 1;


    actions = [3, 4, 3, 4, 1
13, 14, 14, 15, 1
12, 13, 13, 15, 1
11, 12, 12, 15, 1
13, 14, 12, 14, 1
10, 11, 9, 11, 1
8, 15, 6, 15, 1
10, 11, 8, 11, 1
13, 14, 10, 14, 1
11, 12, 9, 12, 1
11, 12, 10, 12, 1
9, 15, 7, 15, 1
11, 12, 9, 12, 1
12, 13, 10, 13, 1
12, 13, 11, 13, 1
10, 15, 8, 15, 1
12, 13, 10, 13, 1
13, 14, 12, 14, 1
11, 15, 9, 15, 1
13, 14, 11, 14, 1
12, 15, 10, 15, 1
13, 15, 11, 15, 1
14, 15, 12, 15, 1
5, 6, 5, 6, 2
14, 15, 15, 16, 1
13, 14, 14, 16, 1
12, 13, 11, 13, 1
10, 16, 8, 16, 1
12, 13, 10, 13, 1
13, 14, 11, 14, 1
13, 14, 12, 14, 1
11, 16, 9, 16, 1
13, 14, 11, 14, 1
14, 15, 13, 15, 1
12, 16, 10, 16, 1
14, 15, 12, 15, 1
13, 16, 11, 16, 1
14, 16, 12, 16, 1
15, 16, 13, 16, 1
7, 8, 7, 8, 2
15, 16, 16, 17, 1
14, 15, 13, 15, 1
12, 17, 10, 17, 1
14, 15, 12, 15, 1
15, 16, 14, 16, 1
13, 17, 11, 17, 1
15, 16, 13, 16, 1
14, 17, 12, 17, 1
15, 17, 13, 17, 1
16, 17, 14, 17, 1
9, 10, 9, 10, 2
16, 17, 15, 17, 1
14, 18, 12, 18, 1
16, 17, 14, 17, 1
15, 18, 13, 18, 1
16, 18, 14, 18, 1
17, 18, 15, 18, 1
11, 12, 11, 12, 2
16, 19, 14, 19, 1
17, 19, 15, 19, 1
18, 19, 16, 19, 1
13, 14, 13, 14, 2
18, 20, 16, 20, 1
19, 20, 17, 20, 1
15, 16, 15, 16, 2
20, 21, 18, 21, 1
17, 18, 17, 18, 1
19, 20, 19, 20, 1
21, 22, 21, 22, 1];


    global bst_cnt;

    % fprintf("bst initiated!!! %d\n", bst_cnt);
    bst_cnt = bst_cnt+1;
    global current_max_depth;
    global depth_sen;
    global group_search_cnt;

    current_max_depth = 0;
    tmp = bst(M, actions, 0, 9, th, 1);

    % param = {M, actions, 0, 9, th};
    % tmp = bst_opt(param);
    if current_max_depth > 0
        % fprintf("bst with depth %d finished.\n", current_max_depth);
        depth_sen = 1;
    end

    M = tmp{1};
    optimizing_param = tmp{2};
    depth = tmp{3};

    global max_depth_in_search;
    if depth > max_depth_in_search
        max_depth_in_search = depth;
    end
    
    if depth > 0
        group_search_cnt = group_search_cnt+1;
    end

    % t2 = toc;
    % keyboard;
    
    

end






