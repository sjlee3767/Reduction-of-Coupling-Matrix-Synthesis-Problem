 % gets a matrix M, pivot (i, j), and the rotation angle theta
 % performs Givens similarity transformation with the given input
 % returns the transformed matrix M' and the rotation matrix R

function [ret, R] = M_rotation(M, i, j, theta)
    R = eye(length(M));
    R(i, i) = cos(theta);
    R(j, j) = cos(theta);
    R(i, j) = -sin(theta);
    R(j, i) = sin(theta);
    ret = R.'*M*R;
end