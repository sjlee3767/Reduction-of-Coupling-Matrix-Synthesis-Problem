function ret = M_trim_zeros(M) % reduces very small elements of the given matrix M to zero
    s = size(M);
    for i = 1:s(1)
        for j = 1:s(2)
            if(abs(real(M(i, j))) < 1e-12)
                M(i, j) = 1j*imag(M(i, j));
            end
            if(abs(imag(M(i, j))) < 1e-12)
                M(i, j) = real(M(i, j));
            end
        end
    end
    ret = M;
end