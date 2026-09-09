function [optimizing_param, M, R] = extended_box_nth_param_4(M_init, th)
    % tic;
    M = M_init;

    N = length(M);
    
    I = eye(N);

    % T = eye(N);
    % T(1, 1) = 0;
    % T(N, N) = 0;
    % T(1, 2) = 1;
    % T(2, 1) = 1;
    % T(N-1, N) = 1;
    % T(N, N-1) = 1;
    % for i = 2:2:N-2
    %     T(i, i+1) = 1;
    %     T(i+1, i) = 1;
    % end
    % for i = 2:N-3
    %     T(i, i+2) = 1;
    %     T(i+2, i) = 1;
    % end
    global target_topology;
    T = target_topology;
    global opt_var_cnt;
    eq_cnt = opt_var_cnt;

    
    P = I;
    
    R = zeros(N);
    R(:, 1) = I(:, 1);
    R(:, 2) = I(:, 2);
    R(:, end-1) = I(:, end-1);
    R(:, end) = I(:, end);
    
    R(:, N-2) = zeros(N, 1);
    % R(4:N-2, N-2) = get_unit_vec_spherical(th);
    % R(4:N-2, N-2) = normalize_C(th);
    global initial_point;
    global tangent;
    tmp  = initial_point' + tangent*th;
    R(4:N-2, N-2) = normalize_C(tmp);

    V = [R(:, 1), R(:, 2), R(:, end-2), R(:, end-1), R(:, end)];
    for i = N-3:-1:3
        R(:, i) = normalize_C(get_next_orthdcp([V], M*R(:, i+2)));
        
        V = [V, R(:, i)];
    end
    M = R'*M*R;

    M = M_trim_zeros(M);

    % tmp = (M.*(ones(N)-T_target)).^2;
    % optimizing_param = sqrt(sum(tmp(:)));


    optimizing_param = zeros(1, eq_cnt);
    for i = 1:N/2-2
        optimizing_param(i) = M(i*2+1, i*2+2)^2;
    end
    for i = 1:N-8
        optimizing_param(i+N/2-2) = M(2, 4+i)^2;
    end
    % for i = 1:N
    %     for j = i:N
    %         if T(i, j) == 0
    %             optimizing_param = [optimizing_param, M(i, j)^2];
    %         end
    %     end
    % end
    % disp(toc);
    % keyboard;
end
