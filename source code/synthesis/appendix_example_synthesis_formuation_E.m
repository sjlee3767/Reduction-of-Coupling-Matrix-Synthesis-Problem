sympref('FloatingPointOutput',true);

% sympref('default');

M = readmatrix("../M.txt");
LR = {}; % put ground truth here for comparison

M_init = M;
global bst_cnt;
bst_cnt=0;

n = length(M)-2;

var_cnt = 2;

% [opts, func, iter_cnt, eps, box_cmp_eps] = set_options(n);

opts = optimoptions('lsqcurvefit');
opts.FunctionTolerance = 1e-10;
opts.OptimalityTolerance = 1e-12;
opts.StepTolerance = 1e-12;
% opts.FunctionTolerance = 1e-20;
% opts.OptimalityTolerance = 1e-22;
% opts.StepTolerance = 1e-22;
opts.MaxIterations = 500;
opts.OutputFcn = @outfunc;
% opts.Algorithm = "trust-region-reflective";
opts.Algorithm = "levenberg-marquardt";
% opts.Display = "iter-detailed";
opts.Display = "off";
iter_cnt = 10000;

func = @appendix_ex_1;

eps = 1e-4;
box_cmp_eps = 0.001;

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

% global idx_type;
% idx_type = 0;

global branch;
tic;
% cnt_tmp = 0;
% sol_cnt2 = zeros(1, 2^(n/2-1));

% sol_cnt_per_idx_type = zeros(2, 2^(n/2-1));
% init_test_cnt = 100;
% init_cnt = 0;
% init_state = 1;

global initial_point;
global tangent;
sol_initial_points = [];

for i = 1:iter_cnt*1000000

    if mod(i, 100) == 1
        elapsed_time = toc; 
        fprintf("%dth iteration done. Found solution (multiplicity) count: %d (%d). Minimum optimizing param found: %f. Miss Count: %d. Elapsed time: %f\n", i, size(sol, 1), sol_cnt_multiplicity, min_opt_param, miss_cnt, elapsed_time);
        % fprintf("%dth iteration done. Found solution (multiplicity) count: %d (%d). Minimum optimizing param found: %f. idx type: %d. Elapsed time: %f\n", i, size(sol, 1), sol_cnt_multiplicity, min_opt_param, idx_type, elapsed_time);
    end

    % th0 = rand(1, var_cnt)*2*pi;
    % th0 = [rand*2*pi rand(1, var_cnt-1)*pi];
    % th0 = [(rand(1, var_cnt-1)-0.5)*pi, rand*2*pi];
    % th0 = randn(1, var_cnt+1);
    % th0 = vec_to_angles(randn(1, var_cnt+1)); % uniform sample of unit vectors

    initial_point = randn(1, var_cnt+1); % uniform sample of unit vectors
    initial_point = normalize_C(initial_point);
    tangent = null(initial_point);
    th0 = zeros(var_cnt, 1);


    % initial_point = adaptive_sampling(var_cnt+1, 100, sol_initial_points);
    % tangent = null(initial_point);
    % th0 = zeros(var_cnt, 1);

    % th0 = [2.095173855071009   4.199084583102414   2.825422594692446];
    % th0 = th0 + (rand(1, var_cnt)-0.5)*0.01;

    % th0 = [   5.591898533472772   3.130409317806429   4.523353333919114   2.023315742151656   0.055646311103160                   0   6.066332592850055   4.771394895398240 5.139902001687362];
    % th0 = [   0.135143246020117   0.150811381194793   3.896754522674184   2.987399776794193   4.367055095652868   3.666208882751067   4.542848424853341   0.878266888084047 5.306256064397186];
    % th0 = [   0.877942397966348   5.534211211920007   2.293231095971694   2.524879068122371   1.799686186836692   0.452015822341190   5.590221108185504   1.570580971065361 2.822471396749281];
    % th0 = [                   0   3.131908012702490   5.490019988249695   6.263638404695647   3.044057600476238   6.257514690883869   4.950600365618175   6.257908741794395 4.712451156278145];
    % th0 = [0.345096290067505   4.554964630819250   4.293249787148342   0.040412953191646   5.844253926089571   3.801705677796647   5.107257417884126   4.203895308753359 4.598152816341628];
    
    target = zeros(1, 2);
    % lb = zeros(1, var_cnt);
    % 
    % rb = ones(1, var_cnt)*2*pi;

    % rb = ones(1, var_cnt)*pi;
    % rb(1) = 2*pi;

    % lb = -ones(1, var_cnt)/2*pi;
    % lb(end) = 0;
    % rb = -lb;
    % rb(end) = 2*pi;

    % lb = -10*ones(1, var_cnt+1);
    % rb = 10*ones(1, var_cnt+1);

    lb = -10*ones(1, var_cnt);
    rb = 10*ones(1, var_cnt);
    % lb = -ones(1, var_cnt)/2;
    % rb = ones(1, var_cnt)/2;
    for k = 1:6
        branch = k;
        depth_sen = 0;
        
        % do numerical optimization
        % th_sol = lsqcurvefit(@(th, M_fixed) func(M_fixed, th), th0, M_init, target, lb, rb, opts);    
        [th_sol, resnorm, residual, exitflag, output] = lsqcurvefit(@(vars, M_fixed) func(M_fixed, vars), th0, M_init, target, lb, rb, opts);   
        exitflagcnt(exitflag+3) = exitflagcnt(exitflag+3)+1;
    
        [optimizing_param, M, R, depth] = func(M_init, th_sol);
    
        tmp  = initial_point' + tangent*th_sol;
        th_sol = normalize_C(tmp)';
    
        % keyboard;
    
        if depth_sen == 1
            fprintf("lsqcurvefit with positive depth finished!\n");
        end
    
        % val = max(optimizing_param.^(1/2), [], "all");
        val = max(abs(optimizing_param), [], "all");
        % pseudo_sol = [pseudo_sol; th_sol];
    
        if val < min_opt_param
            min_opt_param = val;
            min_M = M;
        end
    
        if val < eps % check if actual solution. (global minimum)
            % [M, legit] = fix_boxes(M, box_first, box_cnt, eps);
            % if legit == 0 % check if actual solution. (topology)
            %     sol2 = [sol2; th_sol];
            %     sol_M2{end+1} = M;
            %     miss_cnt = miss_cnt+1;
            %     % keyboard;
            %     continue;
            % end

            % keyboard;
    
            sol_cnt_multiplicity = sol_cnt_multiplicity+1;
            sol_initial_points = [sol_initial_points; initial_point];
            
            [sol, sol_M] = add_solution(sol, th_sol, sol_M, M, box_cmp_eps, LR);
    
            if size(sol, 1) > sol_cnt
                sol_cnt = sol_cnt+1;
                % sol_cnt_per_idx_type(1, idx_type+1) = sol_cnt_per_idx_type(1, idx_type+1)+1;
                % keyboard;
            end
            if depth > 0
                disp("solution with depth > 0 found!");
            end
        end

    end
    
    % idx_type = mod(idx_type+1, 2^(n/2-1));

    if toc >= 30
        break;
    end

end

[sol, sort_idx] = sortrows(sol, var_cnt+2, "descend");
sol_M = sol_M(sort_idx);
sol_M2 = sol_M;

elapsed_time = toc; 
fprintf("Numerical optimization done.\nElapsed time: %f seconds.\nNumber of found solutions: %d.\nNumber of iterations done: %d\n", elapsed_time, size(sol, 1), i);

% for i = 1:size(sol, 1)
%     writematrix(sol_M{i}, sprintf("C:\\sjlee\\Similarity_transformation\\CM reduction test\\data\\extended_box_%dth\\M_%dth_param1_%d", n, n, i));
% end

% for i = 1:length(sol)
%     disp(get_last_couplings(sol_M{i}));
% end

function [sol, sol_M] = add_solution(sol, th, sol_M, M, eps, LR)
    sen = 1;
    if length(sol) == 0
        sol = [th, 1];
        sol_M = {M};
    else
        n = size(sol, 1);
        for i = 1:n
            % tmp = (sol(i,:)-th)';
            % dist = sqrt(tmp'*tmp);
            % if dist < 0.001
            %     sen = 0;
            % end
            if comp_matrices(sol_M{i}, M, eps) == 1
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
        if comp_matrices(LR{i}, M, eps) == 1
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
    % Set absolute threshold on the squared norm of the residual vector
    eps = 1e-4;  % Your desired threshold
    eps = 1e-8;

    % Stop if current residual norm squared is below threshold
    stop = false;
    if strcmp(state, 'iter')
        % disp("@");
        % disp(optimValues.residual');
        % disp(size(optimValues.residual'));
        % disp(optimValues.fval);
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
