sympref('FloatingPointOutput',true);

M = readmatrix("../M.txt");

M_init = M;

N = length(M);
n = N-2;

[opts, func, iter_cnt, eps, box_cmp_eps, T] = set_options(n);

var_cnt = (n-1)*n/2;
eq_cnt = 0;

for i = 1:N
    for j = i:N
        if T(i, j) == 0
            eq_cnt = eq_cnt+1;
        end
    end
end


sol = [];
sol_M = {};
sol_cnt_multiplicity = 0;
sol_cnt = 0;

min_opt_param = 1e10;
min_M = M;

global depth_sen;
global max_depth_in_search;
max_depth_in_search = 0;

tic;
for i = 1:iter_cnt
    if mod(i, 100) == 1
        elapsed_time = toc; 
        fprintf("%dth iteration done. Found solution count: %d. Minimum optimizing param found: %f. Elapsed time: %f\n", i, size(sol, 1), min_opt_param, elapsed_time);
    end

    th0 = rand(1, var_cnt)*2*pi;
    target = zeros(1, eq_cnt);
    lb = zeros(1, var_cnt);
    rb = ones(1, var_cnt)*2*pi;
    
    % do numerical optimization
    th_sol = lsqcurvefit(@(th, M_fixed) func(M_fixed, th), th0, M_init, target, lb, rb, opts);    
    [optimizing_param, M, R] = func(M_init, th_sol);

    val = max(optimizing_param.^(1/2), [], "all");

    if val < min_opt_param
        min_opt_param = val;
        min_M = M;
    end

    
    if val < eps % check if actual solution. (global minimum)
        sol_cnt_multiplicity = sol_cnt_multiplicity+1;
        [sol, sol_M] = add_solution(sol, th_sol, sol_M, M, box_cmp_eps);

        if size(sol, 1) > sol_cnt
            sol_cnt = sol_cnt+1;
        end
    end

    % if toc >= 600
    %     break;
    % end
end

[sol, sort_idx] = sortrows(sol);
sol_M = sol_M(sort_idx);

elapsed_time = toc; 
fprintf("Numerical optimization done.\nElapsed time: %f seconds.\nNumber of found solutions: %d.\n", elapsed_time, size(sol, 1));


function [sol, sol_M] = add_solution(sol, th, sol_M, M, eps)
    if length(sol) == 0
        sol = [th];
        sol_M = {M};
        return;
    end
    n = size(sol, 1);
    sen = 1;
    for i = 1:n
        if comp_matrices(sol_M{i}, M, eps) == 1
            sen = 0;
        end
    end
    if sen == 1
        sol = [sol; th];
        sol_M{end+1} = M;
    else
        
    end
end

function [opts, func, iter_cnt, eps, box_cmp_eps, T] = set_options(n)
    switch n
        case 8
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-12;
            % opts.OptimalityTolerance = 1e-8;
            % opts.StepTolerance = 1e-10;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 10000;
            
            func = @extended_box_8th_param_2;
            
            eps = 1e-4;
            box_cmp_eps = 1e-2;
            
            T = [0     1     0     0     0     0     0     0     0     0
                 1     1     1     1     0     0     0     0     0     0
                 0     1     1     0     1     0     0     0     0     0
                 0     1     0     1     1     1     0     0     0     0
                 0     0     1     1     1     0     1     0     0     0
                 0     0     0     1     0     1     1     1     0     0
                 0     0     0     0     1     1     1     0     1     0
                 0     0     0     0     0     1     0     1     1     0
                 0     0     0     0     0     0     1     1     1     1
                 0     0     0     0     0     0     0     0     1     0];
        case 14
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
            iter_cnt = 20000;

            func = @extended_box_14th_param_2;

            eps = 1e-3;
            box_cmp_eps = 1e-2;

            T = [0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 1     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0
                 0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0
                 0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0
                 0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0
                 0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0
                 0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0
                 0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0
                 0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0
                 0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0
                 0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0
                 0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0
                 0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0
                 0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     0
                 0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0];

        case 22
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-80;
            opts.OptimalityTolerance = 1e-80;
            opts.StepTolerance = 1e-40;
            opts.MaxIterations = 500;
            opts.OutputFcn = @outfunc;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 100000;

            func = @extended_box_22th_param_2;

            eps = 1e-3;
            box_cmp_eps = 1e-2;
            
            T = [0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 1     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     1     0     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     0     1     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     1     1     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0];
        otherwise
            disp("No preset option for the given order of matrix!");
    end
end

function stop = outfunc(x, optimValues, state)
    % Set absolute threshold on the squared norm of the residual vector
    eps = 1e-4;  % Your desired threshold
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