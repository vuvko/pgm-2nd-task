clear;

num_points = 10;
n = 256;
k = 32;
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

ch = 1 + q * log2(q) + (1 - q) * log2(1 - q);
display(['Channel ', num2str(ch)]);
display(['Speed ', num2str(k / n)]);
e = mod(binornd(1, q, [n, 1]), 2);
s = mod(H * e, 2);
for l = 0.1:0.1:1
  disp(l)
  disp('------')
  t = 0;
  for it = 0:5
    tic; ldpc_decoding(s, H, q, 'schedule', 'parallel', 'eps', 1e-4, 'max_iter', 200, 'damping', l, 'display', 0);
    t = t + toc;
  end;
  disp(t/5)
  t = 0;
  for it = 0:5
    tic; ldpc_decoding(s, H, q, 'schedule', 'sequential', 'eps', 1e-4, 'max_iter', 200, 'damping', l, 'display', 0);
    t = t + toc;
  end;
  disp(t/5)
end;
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
