function ret = get_box_absmin(M, i)
    ret = min(abs([M(i, i+2) M(i, i+3) M(i+1, i+2) M(i+1, i+3)]));
    % [v, i] = get_smallest_box_elt([M(i, i+2) M(i, i+3); M(i+1, i+2) M(i+1, i+3)]);
    % ret = v;
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