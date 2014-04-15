clear;

num_points = 10;
n = 50;
k = 10;
m = n - k;
%q = rand(1) / 2; % q --- random
q = 0.1;        % q --- fixed 
H = make_ldpc_mex(m, n, 3);
%t = load('H.mat');
%H = t.H;
% проверяем правильность построения порождающей матрицы mod(H * G, 2) == 0
%[G, ind] = ldpc_gen_matrix(H);
%assert(sum(sum(G(ind, :) ~= eye(k))) == 0);
%assert(sum(sum(mod(H * G, 2) ~= 0)) == 0);
%disp 'ldpc_gen_matrix checked.'

display(['Channel ', num2str(1 + q * log2(q) + (1 - q) * log2(1 - q))]);
display(['Speed ', num2str(k / n)]);
%e = mod(binornd(1, q, [n, 1]), 2);
%s = mod(H * e, 2);
%[e_n, status] = ldpc_decoding(s, H, q, 'schedule', 'parallel', 'eps', 1e-4, 'max_iter', 50, 'damping', 1);
%if status == 2
%    disp 'Max iter.';
%else
%    disp 'Good';
%end
%sum(e ~= e_n) / n

%[err_bit, err_block, diver, b_stat] = ldpc_mc(H, q, num_points);
%plot(b_stat);
%params.max_iter = 200;
%params.damping = 3/8;
%plot_beliefs(H, q, num_points, params);

num_points = 5;
% Шеннон в зависимости от r
n = 256;
q = 0.1;
params.damping = 0.85;
k = [16:16:256];
j = 3;
err_bit = zeros(length(k), 1);
err_block = zeros(length(k), 1);
diver = zeros(length(k), 1);
for i = 1:length(k)
    H = make_ldpc_mex(n - k(i), n, j);
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

figure;
plot(err_bit)
hold on
plot(err_block)
plot(diver)
hold off

% Шеннон в зависимости от n
r = 0.25;
n = [16:16:256];
err_bit = zeros(length(n), 1);
err_block = zeros(length(n), 1);
diver = zeros(length(n), 1);
for i = 1:length(n)
    k = r * n(i);
    H = make_ldpc_mex(n(i) - k, n(i), j);
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

figure;
plot(err_bit)
hold on
plot(err_block)
plot(diver)
hold off

% Шеннон в зависимости от j
n = 256;
k = 32;
m = n - k;
j = [3:1:12];
err_bit = zeros(length(j), 1);
err_block = zeros(length(j), 1);
diver = zeros(length(j), 1);
for i = 1:length(j)
    H = make_ldpc_mex(m, n, j(i));
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

figure;
plot(err_bit)
hold on
plot(err_block)
plot(diver)
hold off
