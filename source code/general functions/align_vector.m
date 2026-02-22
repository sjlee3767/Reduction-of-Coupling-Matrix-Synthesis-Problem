function ret = align_vector(v, w) % puts every elements in v specifed by w into the target axis
    % zero entry in w means that axis should not be rotated.
    n = length(v);
    for i = 1:length(w)
        if w(i) ~= 0 && w(i) ~= 1 && w(i) ~= 2
            disp("topology vector format is not valid");
            return;
        end
        if w(i) == 2
            target_axis = i;
        end
    end
    if length(v) ~= length(w)   % vector lengths must be the same
        disp("dimensions of vectors do not match");
        return;
    end

    P = eye(n);

    % if v is 0-vector in the specified subspace, do nothing.
    sen = 0;
    for i = 1:n
        if w(i) == 0 continue; end
        if v(i) ~= 0 sen = 1; end
    end
    if sen == 0
        ret = eye(n);
        return;
    end

    for i = 1:n
        if w(i) == 0 || i == target_axis   % axis that should not be rotated
            continue; end
        
        rot = eye(n);

        % find the rotation within i and target_axis axes that sets the i-th element of v to 0
        if v(target_axis) == 0
            c = 0;
            s = 1;
        else
            theta = atan(v(i)/v(target_axis));
            c = cos(theta);
            s = sin(theta);
        end
            
        rot(i, i) = c;
        rot(i, target_axis) = -s;
        rot(target_axis, i) = s;
        rot(target_axis, target_axis) = c;

        P = rot*P;
        v = rot*v;

    end

    ret = P;

end