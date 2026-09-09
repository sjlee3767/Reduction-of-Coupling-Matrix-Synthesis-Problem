function [ret, R] = M_rotation(M, i, j, theta)
    % R = eye(length(M));
    % R(i, i) = cos(theta);
    % R(j, j) = cos(theta);
    % R(i, j) = -sin(theta);
    % R(j, i) = sin(theta);
    % ret = R.'*M*R;

    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    tmp = R.'*[M(i, :); M(j, :)];
    % keyboard;
    M(i, :) = tmp(1, :);
    M(j, :) = tmp(2, :);
    tmp = [M(:, i), M(:, j)]*R;
    M(:, i) = tmp(:, 1);
    M(:, j) = tmp(:, 2);
    ret = M;
end