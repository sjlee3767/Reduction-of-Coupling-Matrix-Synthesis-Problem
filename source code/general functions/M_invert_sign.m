% gets a matrix M and the index of the axis to be reversed
% performs similarity transformation with the specified reflection direction
function [ret, R] = M_invert_sign(M, i)
    R = eye(length(M));
    M(i, :) = -M(i, :);
    M(:, i) = -M(:, i);
    R(i, i) = -1;
    ret = M;
end