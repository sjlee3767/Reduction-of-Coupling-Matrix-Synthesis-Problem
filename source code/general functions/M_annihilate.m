% gets a matrix M, pivot (i, j), position to be annihilated (ti, tj), and
% an optional check_multiple variable
% performs a Givens similarity transformation with the specified pivot (i, j),
% that annihilates the specified position (ti, tj).
% when check_multiple == 1 and the given M has a zero-valued to-be-aligned
% vector (under the given pivot and annihilated position), warning message
% is printed.

function [ret, R] = M_annihilate(M, i, j, ti, tj, check_multiple)
    if nargin < 6
        check_multiple = 0;
    end

    sen = 0;
    eps = 1e-6;
    
    if ti == i && tj == i
        theta = atan((-M(i, j) - sqrt(M(i, j)^2 -M(i, i)*M(j, j)))/M(j, j) );
    elseif ti == j && tj == j
        theta = atan((-M(i, j) + sqrt(M(i, j)^2 -M(i, i)*M(j, j)))/M(i, i) );
    elseif ti == i && tj == j
        if abs(M(i, i)-M(j, j)) < eps
            if abs(M(i, j)) < eps
                theta = 0;
                sen = 1;
            else
                theta = pi/4;
            end
        else
            theta = 1/2*atan(2*M(i, j)/(M(i, i)-M(j, j)) );
        end
    elseif ti == i
        if abs(M(j, tj)) < eps && abs(M(i, tj)) < eps
            sen = 1;
        end
        if M(j, tj) == 0
            theta = 0;
        else
            theta = -atan(M(i, tj)/M(j, tj));
        end
    elseif ti == j
        if abs(M(i, tj)) < eps && abs(M(j, tj)) < eps
            sen = 1;
        end
        if M(i, tj) == 0
            % theta = pi/2;
            theta = 0;
        else
            theta = atan(M(j, tj)/M(i, tj));
        end
    elseif tj == i
        if abs(M(ti, j)) < eps && abs(M(ti, i)) < eps
            sen = 1;
        end
        if M(ti, j) == 0
            % theta = pi/2;
            theta = 0;
        else
            theta = -atan(M(ti, i)/M(ti, j));
        end
    elseif tj == j
        if abs(M(ti, i)) < eps && abs(M(ti, j)) < eps
            sen = 1;
        end
        if M(ti, i) == 0
            % theta = pi/2;
            theta = 0;
        else
            theta = atan(M(ti, j)/M(ti, i));
        end
    else
        disp("annihilating position is not within pivot.")
        theta = 0;
    end

    if sen == 1 && check_multiple == 1
        fprintf("one-to-infinite mapping detected! i: %d, j: %d, p: %d, q: %d\n", i, j, ti, tj);
        disp(M_trim_zeros(M));
    end
    
    [ret, R] = M_rotation(M, i, j, theta);
    
end