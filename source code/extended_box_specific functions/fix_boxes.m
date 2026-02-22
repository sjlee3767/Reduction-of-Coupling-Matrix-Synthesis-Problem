% given a modified-shoelace-topology matrix M satisfying the necessary interconnecting couplings condition,
% check if M is actually in extended-box topology, up to change of indices.

function [M, sen] = fix_boxes(M, i, k, eps)
    % i : first box index
    % k : number of boxes

    
    B = M(i:i+1, i+2:i+3);
    if abs(B(2, 1)) < eps || abs(B(2, 2)) < eps
        M = M_rotation(M, i, i+1, pi/2);
    end

    sen = 1;

    i = i+2;
    for cnt = 1:k
        B = M(i-2:i-1, i:i+1);
        if abs(B(1, 2)) < eps
            
        elseif abs(B(1, 1)) < eps
            M = M_rotation(M, i, i+1, pi/2);
        else
            sen = 0;
            break;
        end
        i = i+2;
    end
end

