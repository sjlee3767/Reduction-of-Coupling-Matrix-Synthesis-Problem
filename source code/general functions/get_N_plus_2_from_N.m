function ret = get_N_plus_2_from_N(M, R1, R2)
    n = length(M);
    M = [zeros(n, 1), M, zeros(n, 1)];
    M = [zeros(1, n+2); M; zeros(1, n+2)];
    M(1,2) = sqrt(R1);
    M(2,1) = sqrt(R1);
    M(end,end-1) = sqrt(R2);
    M(end-1,end) = sqrt(R2);
    ret = M;
end