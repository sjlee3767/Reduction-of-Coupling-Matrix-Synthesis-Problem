function ret = get_last_couplings(M)
    ret = [M(2, 3), M(2, 4), M(end-2, end-1), M(end-3, end-1)];
end