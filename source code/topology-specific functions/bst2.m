% same as bst, except for the loss function customized for the appendix
% example.
function ret = bst2(M, actions, depth, th0, valid_sen)
    M_init = M;
    n = size(actions, 1);
    step = 0.01;
    sample_cnt = 30;
    eps = 1e-6;
    eps2 = 1e-4;

    max_depth = depth;
    disp_sen = 0;
    global current_max_depth;

    if current_max_depth < depth
        current_max_depth = depth;
        % fprintf("to-be-aligned is zero. Current max depth: %d\n", current_max_depth);
        disp_sen = 1;
        % keyboard;
    end

    % if depth > 0
    %     disp("depth > 0");
    %     keyboard;
    % end
    
    global bst_cnt;
    global branch;

    for k = 1:n
        arr = actions(k,:);
        i = arr(1); j = arr(2); p = arr(3); q = arr(4); sen = arr(5);
        % sen = 0;
        
        % if sen >= 1 && has_infinite_sol(M, i, j, p, q, eps) == 1
        %     if valid_sen == 0
        %         fprintf("infinite solution, but no hope!\n");
        %     end
        % end

        if valid_sen == 1 && sen >= 1 && has_infinite_sol(M, i, j, p, q, eps) == 1
            mn = 1e10;
            mn_M = M;
            if disp_sen == 1
                % fprintf("to-be-aligned is zero. Current depth: %d, k: %d, bst_cnt: %d, %d %d %d %d\n", depth, k, bst_cnt, i, j, p, q);
                % disp(M_trim_zeros(M));
                % T = array2table(M_trim_zeros(M));
                % disp(T)
            end

            for th = linspace(0, pi, sample_cnt)
                M_tmp = M_rotation(M, i, j, th);

                % keyboard;

                result = bst2(M_tmp, actions(k+1:n,:), depth+1, th0, valid_sen);

                optimizing_param = result{2};
                if optimizing_param'*optimizing_param < mn
                    mn = optimizing_param'*optimizing_param;
                    mn_M = result{1};
                    mn_opt = optimizing_param;
                    max_depth = max(max_depth, result{3});
                end
            end
            M = mn_M;
            ret = {M, mn_opt, max_depth};

            return;
        else
            M = M_annihilate(M, i, j, p, q);
            if sen == 2 % check if the coupling we wish to annihilate through numerical optimization is zero
                if i >= 5
                    if abs(M(i-2, j)) > eps2
                        valid_sen = 0; % the given free variable has no hope of being a solution, so abort.
                    end
                end
            end
            if sen == 3 && branch == 1 % check for the other option for pivot of case 2, which is nonexistent for appendix example 1.
                M = M_rotation(M, i, j, pi/2);
            end
        end
    end

    % optimizing_param = [M(4, 10), M(5, 10)];
    optimizing_param = [M(3, 9), M(4, 9)];

    
    ret = {M, optimizing_param, max_depth};
end


function ret = has_infinite_sol(M, i, j, p, q, eps)
    ret = 0;
    if p == i && q == i
        % not implemented yet
    elseif p == j && q == j
        % not implemented yet
    elseif p == i && q == j
        if abs(M(i, i)-M(j, j)) < eps
            if abs(M(i, j)) < eps
                ret = 1; 
            end
        end
    elseif p == i
        if abs(M(j, q)) < eps && abs(M(i, q)) < eps
            ret = 1;
        end
    elseif p == j
        if abs(M(i, q)) < eps && abs(M(j, q)) < eps
            ret = 1;
        end
    elseif q == i
        if abs(M(p, j)) < eps && abs(M(p, i)) < eps
            ret = 1;
        end
    elseif q == j
        if abs(M(p, i)) < eps && abs(M(p, j)) < eps
            ret = 1;
        end
    else
        disp("annihilating position is not within pivot.")
    end

    
end


function [val, idx] = get_smallest_box_elt(B)
    a = abs(B(1, 1));
    b = abs(B(1, 2));
    c = abs(B(2, 1));
    d = abs(B(2, 2));

    % B = [a b; c d], B = [1 2; 3 4]

    if a <= b && a <= c && a <= d
        val = B(1, 1);
        idx = 1;
    elseif b <= a && b <= c && b <= d
        val = B(1, 2);
        idx = 2;
    elseif c <= a && c <= b && c <= d
        val = B(2, 1);
        idx = 3;
    else
        val = B(2, 2);
        idx = 4;
    end
end