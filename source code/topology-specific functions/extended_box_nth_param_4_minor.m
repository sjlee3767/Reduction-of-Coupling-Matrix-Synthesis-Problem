function [optimizing_param, M, R] = extended_box_nth_param_4_minor(M_init, th)
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
    

    eqn_indices = {};
    for i = 1:N
        tmp = [];
        for j = 1:N
            if T(i, j) == 1
                tmp(end+1) = j;
            end
        end
        eqn_indices{i} = tmp;
    end
    
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

    % V = [R(:, 1), R(:, 2), R(:, end-2), R(:, end-1), R(:, end)];
    for i = N-3:-1:3
        V = [];
        indices = eqn_indices{i+2};
        for j = 2:length(indices)
            V = [V, R(:, indices(j))];
        end
        % if 
        R(:, i) = normalize_C(get_next_orthdcp([V], M*R(:, i+2)));
        
        % V = [V, R(:, i)];
    end
    M = R'*M*R;

    M = M_trim_zeros(M);

    % keyboard;

    % tmp = (M.*(ones(N)-T_target)).^2;
    % optimizing_param = sqrt(sum(tmp(:)));


    optimizing_param = zeros(1, eq_cnt);
    % tmp1 = ((ones(N)-T) .* M(i, j)).^2;
    % tmp2 = (eye(N)-R'*R).^2;
    % optimizing_param = [tmp1(:); tmp2(:)]';

    RtR = R'*R;
    idx = 1;
    for i = 1:N
        for j = i:N
            if T(i, j) == 0
                optimizing_param(idx) = M(i, j)^2;
                idx = idx+1;
            end
        end
    end
    for i = 1:N
        for j = i+1:N
            optimizing_param(idx) = RtR(i, j)^2;
            idx = idx+1;
        end
    end
    % disp(toc);
    
end
