% synthesizes a modified-shoelace topology matrix from the saved
% folded-topolog matrix (M.txt), with random free variable values.

M = readmatrix("../M.txt");

n = length(M)-2;

th = rand(1, n/2-2) * 2*pi;

M = folded_to_sword(M, th);

T1 = array2table(M_trim_zeros(M));

M = sword_to_shoelace(M);

T2 = array2table(M_trim_zeros(M));
T2