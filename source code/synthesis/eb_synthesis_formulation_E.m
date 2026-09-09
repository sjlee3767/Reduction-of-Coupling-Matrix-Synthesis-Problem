sympref('FloatingPointOutput',true);

M = readmatrix("../M.txt");
LR = {}; % put ground truth here for comparison

M_init = M;
global bst_cnt;
bst_cnt=0;

n = length(M)-2;

box_first = 3;
box_cnt = n/2-2;
var_cnt = n/2-2;

[opts, func, iter_cnt, eps, box_cmp_eps] = set_options(n);
global var_type;
var_type = 0; % 0, 1, 2 or 3. 0 corresponds to extended box, and other values for variants
if var_type ~= 0 && n == 14 % variants implementation done for n == 14 only
    func = @extended_box_14th_vars;
end

global pseudo_sol;
sol = [];
sol_M = {};
sol_cnt_multiplicity = 0;
sol_cnt = 0;


min_opt_param = 1e10;
min_M = M;

sol2 = [];
sol_M2 = {};
miss_cnt  = 0;

exitflagcnt = zeros(1, 7);

global depth_sen;
global max_depth_in_search;
max_depth_in_search = 0;

global sol_pair;
sol_pair = [];
for i = 1:length(LR) % change index to match topology
    LR{i} = M_trim_zeros(LR{i});
    for j = 4:4:n-1
        LR{i} = M_rotation(LR{i}, j, j+1, pi/2);
    end
    LR{i} = M_trim_zeros(LR{i});
end

LR2 = {};
for i = 1:length(LR)
    sen = 1;
    for j = i+1:length(LR)
        if comp_matrices(LR{i}, LR{j}, 1e-4) == 1
            LR2{end+1} = LR{i};
        end
    end
end
LR = LR2;

global branch;
tic;
global initial_point;
global tangent;
sol_initial_points = [];

for i = 1:iter_cnt

    if mod(i, 100) == 1
        elapsed_time = toc; 
        fprintf("%dth iteration done. Found solution (multiplicity) count: %d (%d). Minimum optimizing param found: %f. Miss Count: %d. Elapsed time: %f\n", i, size(sol, 1), sol_cnt_multiplicity, min_opt_param, miss_cnt, elapsed_time);
    end

    initial_point = randn(1, var_cnt+1); % uniform sample of unit vectors
    initial_point = normalize_C(initial_point);
    tangent = null(initial_point);
    th0 = zeros(var_cnt, 1);

    target = zeros(1, box_cnt);
    lb = -10*ones(1, var_cnt);
    rb = 10*ones(1, var_cnt);
    for k = 0:1
        branch = k;
        depth_sen = 0;
        
        % do numerical optimization
        [th_sol, resnorm, residual, exitflag, output] = lsqcurvefit(@(vars, M_fixed) func(M_fixed, vars), th0, M_init, target, lb, rb, opts);   
        exitflagcnt(exitflag+3) = exitflagcnt(exitflag+3)+1;
    
        [optimizing_param, M, R, depth] = func(M_init, th_sol);
    
        tmp  = initial_point' + tangent*th_sol;
        th_sol = normalize_C(tmp)';
    
        % keyboard;
    
        if depth_sen == 1
            fprintf("lsqcurvefit with positive depth finished!\n");
        end
    
        val = max(abs(optimizing_param).^(1/2), [], "all");
    
        if val < min_opt_param
            min_opt_param = val;
            min_M = M;
        end
    
        if val < eps % check if actual solution. (global minimum)
            sol_cnt_multiplicity = sol_cnt_multiplicity+1;
            sol_initial_points = [sol_initial_points; initial_point];
            
            [sol, sol_M] = add_solution(sol, th_sol, sol_M, M, box_first, box_cnt, box_cmp_eps, LR);
    
            if size(sol, 1) > sol_cnt
                sol_cnt = sol_cnt+1;
            end
            if depth > 0
                disp("solution with depth > 0 found!");
            end
        end

    end

    if toc >= 600
        break;
    end

    % if sol_cnt == 8
    %     break;
    % end
       

end

[sol, sort_idx] = sortrows(sol, var_cnt+2, "descend");
sol_M = sol_M(sort_idx);
sol_M2 = sol_M;

elapsed_time = toc; 
fprintf("Numerical optimization done.\nElapsed time: %f seconds.\nNumber of found solutions: %d.\nNumber of iterations done: %d\n", elapsed_time, size(sol, 1), i);

function [sol, sol_M] = add_solution(sol, th, sol_M, M, box_first, box_cnt, eps, LR)
    sen = 1;
    if length(sol) == 0
        sol = [th, 1];
        sol_M = {M};
    else
        n = size(sol, 1);
        for i = 1:n
            if comp_boxes(sol_M{i}, M, box_first, box_cnt, eps) == 1
                sol(i, end) = sol(i, end)+1;
                sen = 0;
            end
        end
        if sen == 1
            sol = [sol; th, 1];
            sol_M{end+1} = M;
        else
            % disp("!");
        end
    end

    m = length(LR);
    sen2 = 1;
    idx = -1;
    for i = 1:m
        if comp_boxes(LR{i}, M, box_first, box_cnt, eps) == 1
            sen2 = 0;
            idx = i;
        end
    end
    if sen == 1 && sen2 == 1
        disp("new solution found, but not in LR!");
        % keyboard;
    elseif sen == 1 && sen2 == 0
        global sol_pair;
        sol_pair = [sol_pair, idx];
        disp("new valid solution found!");
    end
end

function stop = outfunc(x, optimValues, state)
    eps = 1e-4;
    eps = 1e-8;

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

function [opts, func, iter_cnt, eps, box_cmp_eps] = set_options(n)
    switch n
        case 8
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-12;
            % opts.OptimalityTolerance = 1e-12;
            % opts.StepTolerance = 1e-14;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 1000;

            func = @extended_box_8th;

            eps = 1e-4;
            box_cmp_eps = 0.001;
        case 10
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-12;
            opts.OptimalityTolerance = 1e-12;
            opts.StepTolerance = 1e-14;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 1000000;

            func = @extended_box_10th_test;

            eps = 1e-4;
            box_cmp_eps = 0.001;

        case 14
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-12;
            opts.OptimalityTolerance = 1e-8;
            opts.StepTolerance = 1e-14;
            % opts.FunctionTolerance = 1e-10;
            % opts.OptimalityTolerance = 1e-8;
            % opts.StepTolerance = 1e-12;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 2000000;
            
            func = @extended_box_14th;
            
            eps = 1e-3;
            box_cmp_eps = 0.001;
            
        case 22
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-40;
            opts.OptimalityTolerance = 1e-40;
            opts.StepTolerance = 1e-80;
            % opts.FunctionTolerance = 1e-18;
            % opts.OptimalityTolerance = 1e-18;
            % opts.StepTolerance = 1e-20;
            % opts.FunctionTolerance = 1e-12;
            % opts.OptimalityTolerance = 1e-8;
            % opts.StepTolerance = 1e-10;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 100000;

            func = @extended_box_22th;

            eps = 1e-3;
            box_cmp_eps = 0.01;
        otherwise
            disp("No preset option for the given order of matrix!");
    end
end