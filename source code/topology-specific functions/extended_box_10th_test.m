function [optimizing_param, M, R, depth] = extended_box_10th_test(M_init, th)
    global initial_point;
    global tangent;

    M = M_init;

    N = length(M);
    n = N-2;

    if length(th) ~= n/2-2
        fprintf("number of free variables does not match the order of M!\n");
    end

    I = eye(N);
    P = I;

    r0 = I(:, 1);
    r1 = I(:, 2);
    r2 = I(:, 3);
    r9 = I(:, end-2);
    r10 = I(:, end-1);
    r11 = I(:, end);

    r8 = normalize_C(get_next_orthdcp([r9 r10], M*r9));

    P1 = align_vector(M*r8, [0 0 0 1 1 1 1 2 0 0 0 0]);
    P = P1*P;

    Pr3 = zeros(N, 1);
    tmp = initial_point' + tangent*th;
    Pr3(4:7) = normalize_C(tmp);
    r3 = P'*Pr3;

    v1 = normalize_C((I-r1*r1'-r2*r2'-r3*r3')*M*r2);
    v2 = normalize_C((I-r2*r2'-r3*r3')*M*r3);

    v3 = normalize_C((I-v2*v2')*v1); % {v2, v3} is an orthogonal basis of the linear space of {v1, v2}
    Mr8 = M*r8;
    r5 = normalize_C((Mr8'*v3)*v3 + (Mr8'*v2)*v2);
    r4 = normalize_C(get_next_orthdcp([r1 r2 r3 r5], M*r2));

    % P2 = align_vector(P*r3, [0 0 0 2 1 1 1 1 0 0 0 0]);
    % P = P2*P;
    % P3 = align_vector(P*v1, [0 0 0 0 2 1 1 1 0 0 0 0]);
    % P = P3*P;
    % P4 = align_vector(P*v2, [0 0 0 0 0 2 1 1 0 0 0 0]);
    % P = P4*P;
    % P5 = align_vector(P*M*r8, [0 0 0 0 1 2 0 0 0 0 0 0]);
    % P = P5*P;
    % Pr4 = [0 0 0 0 1 0 0 0 0 0 0 0]';
    % Pr5 = [0 0 0 0 0 1 0 0 0 0 0 0]';
    % r4 = P'*Pr4;
    % r5 = P'*Pr5;

    r6 = normalize_C(get_next_orthdcp([r2 r3 r4 r5], M*r4));
    r7 = normalize_C(get_next_orthdcp([r2 r3 r4 r5 r6 r8], M*r5));

    R = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11];
    M = R'*M*R;

    actions = [8, 9, 6, 9, 1
               3, 4, 3, 4, 3
               5, 6, 3, 6, 1
               7, 8, 5, 8, 1
               9, 10, 7, 10, 1];

    tmp = bst(M, actions, 0, n/2-2, th, 1);
    M = tmp{1};
    optimizing_param = tmp{2};
    depth = tmp{3};

    M = M_trim_zeros(M);
end
