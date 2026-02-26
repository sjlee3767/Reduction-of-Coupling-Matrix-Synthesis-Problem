sympref('FloatingPointOutput',true);

M = readmatrix("../M.txt");
M_init = M;

global bst_cnt;
bst_cnt=0;

% 8th
opts = optimoptions('lsqcurvefit');
opts.FunctionTolerance = 1e-12;
opts.MaxIterations = 500;
opts.OutputFcn = @outfunc;
opts.Algorithm = "levenberg-marquardt";
opts.Display = "off";
iter_cnt = 1000;

func = @extended_box_8th_C;

box_first = 3;
box_cnt = 2;
var_cnt = 2;

eps = 1e-4;
box_cmp_eps = 0.001;

% 14th
% opts = optimoptions('lsqcurvefit');
% opts.FunctionTolerance = 1e-12;
% opts.OptimalityTolerance = 1e-12;
% opts.StepTolerance = 1e-14;
% opts.MaxIterations = 500;
% opts.OutputFcn = @outfunc;
% opts.Algorithm = "levenberg-marquardt";
% opts.Display = "off";
% iter_cnt = 10000;
% 
% func = @extended_box_14th_C;
% 
% box_first = 3;
% box_cnt = 5;
% var_cnt = 5;
% 
% eps = 1e-4;
% box_cmp_eps = 0.001;

% 22nd
% opts = optimoptions('lsqcurvefit');
% opts.FunctionTolerance = 1e-80;
% opts.OptimalityTolerance = 1e-80;
% opts.StepTolerance = 1e-40;
% opts.MaxIterations = 500;
% opts.OutputFcn = @outfunc;
% opts.Algorithm = "levenberg-marquardt";
% opts.Display = "off";
% iter_cnt = 100000;
% 
% func = @extended_box_22th_C;
% 
% box_first = 3;
% box_cnt = 9;
% var_cnt = 9;
% 
% eps = 1e-3;
% box_cmp_eps = 0.01;


sol = [];
sol_M = {};
sol_cnt_multiplicity = 0;
sol_cnt = 0;

min_opt_param = 1e10;
min_M = M;

sol2 = [];
sol_M2 = {};
miss_cnt  = 0;

global depth_sen;
global max_depth_in_search;
max_depth_in_search = 0;

tic;
for i = 1:iter_cnt
    if mod(i, 2) == 1
        elapsed_time = toc; 
        fprintf("%dth iteration done. Found solution (multiplicity) count: %d (%d). Minimum optimizing param found: %f. Miss Count: %d. Elapsed time: %f\n", i, size(sol, 1), sol_cnt_multiplicity, min_opt_param, miss_cnt, elapsed_time);
    end

    th0 = rand(1, var_cnt)*2*pi;
    target = zeros(1, box_cnt);
    lb = zeros(1, var_cnt);
    rb = ones(1, var_cnt)*2*pi;

    depth_sen = 0;
    
    % do numerical optimization
    th_sol = lsqcurvefit(@(th, M_fixed) func(M_fixed, th), th0, M_init, target, lb, rb, opts);    
    [optimizing_param, M, R, depth] = func(M_init, th_sol);

    if depth_sen == 1
        fprintf("lsqcurvefit with positive depth finished!\n");
    end

    % val = max(optimizing_param.^(1/2), [], "all");
    val = max(abs(optimizing_param).^(1/2), [], "all");

    if val < min_opt_param
        min_opt_param = val;
        min_M = M;
    end

    if val < eps % check if actual solution. (global minimum)
        [M, legit] = fix_boxes(M, box_first, box_cnt, eps);
        if legit == 0 % check if actual solution. (topology)
            sol2 = [sol2; th_sol];
            sol_M2{end+1} = M;
            miss_cnt = miss_cnt+1;
            continue;
        end

        sol_cnt_multiplicity = sol_cnt_multiplicity+1;
        [sol, sol_M] = add_solution(sol, th_sol, sol_M, M, box_first, box_cnt, box_cmp_eps);
        if depth > 0
            disp("solution with depth > 0 found!");
        end
    end

end

[sol, sort_idx] = sortrows(sol);
sol_M = sol_M(sort_idx);

elapsed_time = toc; 
fprintf("Numerical optimization done.\nElapsed time: %f seconds.\nNumber of found solutions: %d.\n", elapsed_time, size(sol, 1));

function [sol, sol_M] = add_solution(sol, th, sol_M, M, box_first, box_cnt, eps)
    if length(sol) == 0
        sol = [th];
        sol_M = {M};
        return;
    end
    n = size(sol, 1);
    sen = 1;
    for i = 1:n
        if comp_boxes(sol_M{i}, M, box_first, box_cnt, eps) == 1
            sen = 0;
        end
    end
    if sen == 1
        sol = [sol; th];
        sol_M{end+1} = M;
    else
        % disp("!");
    end
end

function stop = outfunc(x, optimValues, state)
    % Set absolute threshold on the squared norm of the residual vector
    eps = 1e-4;  % Your desired threshold

    % Stop if current residual norm squared is below threshold
    stop = false;
    if strcmp(state, 'iter')
        val = max(abs(optimValues.residual).^(1/2), [], "all");
        if val < eps
            stop = true;
        end
        if optimValues.lambda == 0
            disp("lambda is zero. skipping...");
            stop = true;
        end
    end

end
