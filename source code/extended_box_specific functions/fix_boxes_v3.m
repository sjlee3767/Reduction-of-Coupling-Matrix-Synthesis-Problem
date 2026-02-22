% given a modified-shoelace-topology matrix M satisfying the necessary interconnecting couplings condition,
% check if M is actually in extended-box-variant topology (Fig. 8 (c)), up to change of indices.

function [M, sen] = fix_boxes_v3(M, i, k, eps)
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
        if i == 7 || i == 11
            if abs(B(2, 1)) < eps
                
            elseif abs(B(2, 2)) < eps
                M = M_rotation(M, i, i+1, pi/2);
            else
                sen = 0;
                break;
            end
        else
            if abs(B(1, 2)) < eps
                
            elseif abs(B(1, 1)) < eps
                M = M_rotation(M, i, i+1, pi/2);
            else
                sen = 0;
                break;
            end
        end
        i = i+2;
    end
end

