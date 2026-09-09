function [M] = generalized_box_nth(M_init, th)
    % tic;
    M = M_init;

    n = length(M)-2;

    if length(th) ~= n/2-2
        fprintf("number of free variables does not match the order of M!\n");
    end
    
    M = folded_to_sword(M, th);
    M = sword_to_generalized_box(M);
end






