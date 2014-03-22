clear;

n = 256;
k = 128;
m = n - k;
%q = rand(1) / 2; % q --- random
q = 0.01;        % q --- fixed 
H = make_ldpc_mex(m, n, 4);
num_points = 20;
% проверяем правильность построения порождающей матрицы mod(H * G, 2) == 0
%[G, ind] = ldpc_gen_matrix(H);
%assert(sum(sum(G(ind, :) ~= eye(k))) == 0);
%assert(sum(sum(mod(H * G, 2) ~= 0)) == 0);
%disp 'ldpc_gen_matrix checked.'

display(['Channel ', num2str(1 + q * log2(q) + (1 - q) * log2(1 - q))]);
display(['Speed ', num2str(k / n)]);
e = mod(binornd([1:n]', q), 2);
%v = randi(2, n, 1) - 1;
v = ones(n, 1);
w = xor(v, e);
s = mod(H * w, 2);
[e_n, status] = ldpc_decoding(s, H, q, 'schedule', 'sequential', 'eps', 1e-4);
if status == 2
    disp 'Max iter.';
else
    disp 'Good';
    sum(e ~= e_n)
    max(e_n)
end

[err_bit, err_block, diver] = ldpc_mc(H, q, num_points)
