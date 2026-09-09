global filename_prefix;
filename_prefix = "C:\sjlee\images\sword_to_generalized_box_transformations\test\sword2generalizedbox";
M = readmatrix("../M.txt");

sympref('FloatingPointOutput',true);

% sympref('default');

M_init = M;

N = length(M);
n = N-2;

th = rand(1, n/2-2) * 2*pi;
M = folded_to_sword(M, th);
M_sword = M;

% % M =[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1
% % 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0
% % ];
% 
% 
% N = length(M);
% n = N-2;
% M = rand(N, N).*M;
% for i = 1:N
%     for j = i+1:N
%         M(j, i) = M(i, j);
%     end
% end


global cnt;
cnt = 0;

I = eye(N);

l = 7;
r = n/2+3;
k = (r-l)-5;

write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Sword Topology");

while r < n-1
    M = level_1(M, l+1, r+1, k);
    

    l = l + 2;
    r = r + 1;
    k = k - 1;
end

% ret = M;
% return;
M = M_annihilate(M, 3, 4, 3, 4);
write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Last alignment");
for i = 5:2:N-3
    M = M_annihilate(M, i, i+1, i-2, i+1);
    write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Last alignment");
end




function ret = level_1(M, l, r, k) % k+1: number of level_2 to be performed
    global filename_prefix;
    global cnt;

    for i = 0:k
        M = level_2(M, r, r-i-1);
    end
    
    for i = 1:r-l % r-l: number of column elements to be annihilated
        if (r)-(l+i-1) > 3
            M = M_annihilate(M, l+1+i, l+2+i, l+i, l+2+i);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 1");
        end
        
        
        M = M_annihilate(M, l+i-1, r, l-3+i, r);
        write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 1");
        
        if (r)-(l+i-1) > 3
            M = M_annihilate(M, l+1+i, l+2+i, l-1+i, l+2+i);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 1");
        end
        
        if i < k+1 % do level_2_inv only when level_2 was actually performed
            M = level_2_inv(M, r, (r-k-1)+i-1);
        end
    end

    ret = M;

end

function ret = level_2(M, r, k) % k: index of row to be annihilated
    global filename_prefix;
    global cnt;

    M = M_annihilate(M, k-1, k, k, r);
    write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2");
    for i = 1:r-k-2
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, k, r-i);
        write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2");
        M = level_3_inv(M, r, i);
    end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, r-i-2, r-i);
        write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2");
        M = level_3_inv(M, r, i);
     end

     ret = M;
end

function ret = level_2_inv(M, r, k) % k: index of row to be revived
    global filename_prefix;
    global cnt;

     for i = 1:r-k-2
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, k-2, r-i);
        write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2 Rev");
        M = level_3_inv(M, r, i); 
     end

     for i = r-k-3:-1:1
        M = level_3(M, r, i);
        M = M_annihilate(M, r-i-1, r-i, r-i-2, r-i);
        write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2 Rev");
        M = level_3_inv(M, r, i);
     end

     M = M_annihilate(M, k-1, k, k-3, k);
     write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 2 Rev");

     ret = M;
end

function ret = level_3(M, r, k)
    global filename_prefix;
    global cnt;

    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = 3:2:k
            M = M_annihilate(M, r-i+1, r-i+2, r-i+1, r-i+2);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 3");
        end
    else
        for i = 4:2:k
            M = M_annihilate(M, r-i+1, r-i+2, r-i+1, r-i+2);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 3");
        end
    end
    ret = M;
end

function ret = level_3_inv(M, r, k)
    global filename_prefix;
    global cnt;

    if k <= 2
        ret = M;
        return;
    end
    if mod(k, 2) == 1
        for i = k:-2:3
            M = M_annihilate(M, r-i+1, r-i+2, r-i-1, r-i+2);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 3 Rev");
        end
    else
        for i = k:-2:4
            M = M_annihilate(M, r-i+1, r-i+2, r-i-1, r-i+2);
            write_topology(M, filename_prefix+sprintf("%03d", cnt)+".txt", "Level 3 Rev");
        end
    end
    ret = M;
end

function write_topology(M, filename, caption)
    global cnt;
    T = M;
    n = length(T);
    for i = 1:n
        for j = 1:n
            if abs(T(i, j)) < 1e-12
                T(i, j) = 0;
            else
                T(i, j) = 1;
            end
        end
    end
    fp = fopen(filename, "w");
    fprintf(fp, "%d\n%s\n", cnt, caption);
    fclose(fp);
    writematrix(T, filename, 'WriteMode', 'append');
    cnt = cnt+1;
end