M = readmatrix("../M.txt");

sympref('FloatingPointOutput',true);

% sympref('default');

M_init = M;

N = length(M);
n = N-2;

I = eye(N);

l = 7;
r = n/2+3;
k = (r-l)-5;

p = 3;
q = 1;
while r < n-1
    if p >= 5
        M = output(M, p, p+1, p, p+1, 2);
    else
        M = output(M, p, p+1, p, p+1);
    end
    M = level_1(M, l+1, r+1, k);

    if p == 5
        [v, idx] = get_smallest_box_elt(M(p-2:p-1, p:p+1));
        if idx == 3 || idx == 4
            % M = M_rotation(M, p-2, p-1, pi/2);
        end
    end
    if p >= 5
        [v, idx] = get_smallest_box_elt(M(p-2:p-1, p:p+1));
        % if idx == 1
        %     M = M_rotation(M, p, p+1, pi/2);
        % end
        if abs(M(p-2, p)) < abs(M(p-2, p+1))
            % M = M_rotation(M, p, p+1, pi/2);
        end
        optimizing_param(q) = v;
        q = q+1;
    end

    T = array2table(M_trim_zeros(M));
    % keyboard;

    

    l = l + 2;
    r = r + 1;
    k = k - 1;
    p = p+2;
end

% ret = M;
% return;
for i = p:2:N-3
    M = output(M, i, i+1, i, i+1);
    [v, idx] = get_smallest_box_elt(M(i-2:i-1, i:i+1));
    if idx == 1
        % M = M_rotation(M, i, i+1, pi/2);
    end
    optimizing_param(q) = v;
    q = q+1;
end

T = array2table(M_trim_zeros(M));
% keyboard;

depth = 1;

function ret = level_1(M, l, r, k) % k+1: number of level_2 to be performed
    for i = 0:k
        M = level_2(M, r, r-i-1);
    end
    
    for i = 1:r-l % r-l: number of column elements to be annihilated
        if (r)-(l+i-1) > 3
            M = output(M, l+1+i, l+2+i, l+i, l+2+i);
        end
        
        
        M = output(M, l+i-1, r, l-3+i, r);
        
        if (r)-(l+i-1) > 3
            M = output(M, l+1+i, l+2+i, l-1+i, l+2+i);
        end
        
        if i < k+1 % do level_2_inv only when level_2 was actually performed
            M = level_2_inv(M, r, (r-k-1)+i-1);
        end
    end

    ret = M;

end

function ret = level_2(M, r, k) % k: index of row to be annihilated
    M = output(M, k-1, k, k, r);
    for i = 1:r-k-2
        M = level_3(M, r, i);
        M = output(M, r-i-1, r-i, k, r-i);
        M = level_3_inv(M, r, i);
    end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = output(M, r-i-1, r-i, r-i-2, r-i);
        M = level_3_inv(M, r, i);
     end

     ret = M;
end

function ret = level_2_inv(M, r, k) % k: index of row to be revived
     for i = 1:r-k-2
        M = level_3(M, r, i);
        M = output(M, r-i-1, r-i, k-2, r-i);
        M = level_3_inv(M, r, i); 
     end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = output(M, r-i-1, r-i, r-i-2, r-i);
        M = level_3_inv(M, r, i);
     end

     M = output(M, k-1, k, k-3, k);

     ret = M;
end

function ret = level_3(M, r, k)
    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = 3:2:k
            M = output(M, r-i+1, r-i+2, r-i+1, r-i+2);
        end
    else
        for i = 4:2:k
            M = output(M, r-i+1, r-i+2, r-i+1, r-i+2);
        end
    end
    ret = M;
end

function ret = level_3_inv(M, r, k)
    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = k:-2:3
            M = output(M, r-i+1, r-i+2, r-i-1, r-i+2);
        end
    else
        for i = k:-2:4
            M = output(M, r-i+1, r-i+2, r-i-1, r-i+2);
        end
    end
    ret = M;
end

function [val, idx] = get_smallest_box_elt(B)
    a = abs(B(1, 1));
    b = abs(B(1, 2));
    c = abs(B(2, 1));
    d = abs(B(2, 2));

    % B = [a b; c d], B = [1 2; 3 4]

    if a <= b && a <= c && a <= d
        val = B(1, 1);
        idx = 1;
    elseif b <= a && b <= c && b <= d
        val = B(1, 2);
        idx = 2;
    elseif c <= a && c <= b && c <= d
        val = B(2, 1);
        idx = 3;
    else
        val = B(2, 2);
        idx = 4;
    end
end

function ret =  output(M, i, j, p, q, sen)
    if nargin < 6
        sen = 1;
    end
    fprintf("%d, %d, %d, %d, %d\n", i, j, p, q, sen);
    ret = M;
end