clear;

num_points = 4;
n = 16;
k = 4;
m = n - k;
j = 3;
%q = rand(1) / 2; % q --- random
q = 0.1;        % q --- fixed 
H = make_ldpc_mex(m, n, j);

display(['Channel ', num2str(1 + q * log2(q) + (1 - q) * log2(1 - q))]);
display(['Speed ', num2str(k / n)]);
