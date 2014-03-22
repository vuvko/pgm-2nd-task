clear;

n = 256;
m = 128;
k = n - m;
%q = rand(1) / 2; % q --- random
q = 0.25;        % q --- fixed 
H = make_ldpc_mex(m, n, 3);
% проверяем правильность построения порождающей матрицы mod(H * G, 2) == 0
[G, ind] = ldpc_gen_matrix(H);
assert(sum(sum(G(ind, :) ~= eye(k))) == 0);
assert(sum(sum(mod(H * G, 2) ~= 0)) == 0);
disp 'ldpc_gen_matrix checked.'

e = binornd([1:n]', q);
v = rand(n, 1);
w = v + e;
