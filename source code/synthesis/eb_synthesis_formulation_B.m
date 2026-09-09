sympref('FloatingPointOutput',true);
M = readmatrix("../M.txt");
M_init = M;


N = length(M);
n = N-2;

[opts, func, iter_cnt, eps, box_cmp_eps, T] = set_options(n);

var_cnt = 0;
for i = 1:N
    for j = i:N
        if T(i, j) == 1
            var_cnt = var_cnt+1;
        end
    end
end

eq_cnt = var_cnt; % # of nonzeros of folded = extended_box.

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

    vars0 = (rand(1, var_cnt)-1/2)*2;        
    target = zeros(1, eq_cnt);
    
    % do numerical optimization
    vars_sol = lsqcurvefit(@(vars, M_fixed) func(M_fixed, vars), vars0, M_init, target, [], [], opts);    
    [optimizing_param, M_folded, M] = func(M_init, vars_sol);

    val = max(optimizing_param.^(1/2), [], "all");

    if val < min_opt_param
        min_opt_param = val;
        min_M = M;
    end

    if val < eps % check if actual solution. (global minimum)
        sol_cnt_multiplicity = sol_cnt_multiplicity+1;
        [sol, sol_M] = add_solution(sol, vars_sol, sol_M, M, box_cmp_eps);

        if size(sol, 1) > sol_cnt
            sol_cnt = sol_cnt+1;
        end
    end

    if toc >= 3600*2
        break;
    end

    % if sol_cnt == 8
    %     break;
    % end

end

[sol, sort_idx] = sortrows(sol);
sol_M = sol_M(sort_idx);

elapsed_time = toc; 
fprintf("Numerical optimization done.\nElapsed time: %f seconds.\nNumber of found solutions: %d.\n", elapsed_time, size(sol, 1));

for i = 1:size(sol, 1)
    writematrix(sol_M{i}, sprintf("C:\\sjlee\\Similarity_transformation\\CM reduction test\\data\\extended_box_8th\\M_8th_param3_%d", i));
end


function [sol, sol_M] = add_solution(sol, vars, sol_M, M, eps)
    if length(sol) == 0
        sol = [vars];
        sol_M = {M};
        return;
    end
    n = size(sol, 1);
    sen = 1;
    for i = 1:n
        tmp = (sol(i,:)-vars)';
        dist = sqrt(tmp'*tmp);
        if dist < 0.01
            sen = 0;
        end
        if comp_matrices(sol_M{i}, M, eps) == 1
            sen = 0;
        end
    end
    if sen == 1
        sol = [sol; vars];
        sol_M{end+1} = M;
    else
        
    end
end


function [opts, func, iter_cnt, eps, box_cmp_eps, T] = set_options(n)
    switch n
        case 8
            opts = optimoptions('lsqcurvefit');
            opts.FunctionTolerance = 1e-12;
            % opts.OptimalityTolerance = 1e-12;
            % opts.StepTolerance = 1e-14;
            opts.MaxIterations = 500;
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 1000;
            
            func = @extended_box_8th_param_3;
            
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
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 200000;

            func = @extended_box_14th_param_3;

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
            % opts.Algorithm = "trust-region-reflective";
            opts.Algorithm = "levenberg-marquardt";
            % opts.Display = "iter-detailed";
            opts.Display = "off";
            iter_cnt = 100000;

            func = @extended_box_22th_param_3;

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