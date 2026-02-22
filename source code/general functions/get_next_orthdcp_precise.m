function [ret, V3] = get_next_orthdcp_precise(V1, V2, v) % repeat orthdcp until innner product between the result and all vectors in V are less than criterion
    V = [V1 V2];

    [n, m] = size(V);
    
    ret = normalize_C3(normalize_C3(get_next_orthdcp(V, v)));
    for cnt = 1:400
        sen = 0;
        for i = 1:m
            if abs(ret'*V(:, i)) > 1e-12
                sen = 1;
            end
        end
        if sen == 0
            V3 = [V1, ret];
            return;
        end
        ret = normalize_C3(normalize_C3(get_next_orthdcp(V, ret)));
    end
    V3 = [V1, ret];
    fprintf("warning! get_next_orthdcp_precise could not sufficiently reduce numerical error within given iteration!\n");
end
