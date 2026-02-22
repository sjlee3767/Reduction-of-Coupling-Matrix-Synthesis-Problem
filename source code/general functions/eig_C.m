function [V, LAM] = eig_C(M) % eigendecomposition, with sorted eigenvalues
    [V_tmp, LAM_tmp] = eig(M);
    lam = diag(LAM_tmp);
    [lam_sorted, idx] = sort(lam);
    V = V_tmp;
    LAM = eye(length(V));
    for i = 1:length(V)
        V(:, i) = V_tmp(:, idx(i));
        LAM(i, i) = lam_sorted(i);
    end
end