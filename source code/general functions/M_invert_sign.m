function [ret, R] = M_invert_sign(M, i)
    R = eye(length(M));
    M(i, :) = -M(i, :);
    M(:, i) = -M(:, i);
    R(i, i) = -1;
    ret = M;
end