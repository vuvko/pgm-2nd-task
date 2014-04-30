% строим все нужные графики

init;

h_beliefs_p = figure;
set(h_beliefs_p, 'Color', 'w');
h_beliefs_s = figure;
set(h_beliefs_s, 'Color', 'w');
%damping = 0.1:0.1:1;
%color = [[150, 75, 0]; [235, 0, 0]; [230, 230, 5]; [70, 140, 250]; ...
%    [0, 0, 235]; [220, 115, 150]; [0, 235, 0]; [0, 135, 235]; ...
%    [230, 50, 255]; [50, 240, 50]] / 255.0;
%for i = 1:length(damping)
for damping = 0.1:0.1:1
    params.damping = damping;
    params.schedule = 'parallel';
    params.style = '-';
    params.color = [damping, 1 - damping, 0.0];
    figure(h_beliefs_p);
    plot_beliefs(H, q, num_points, params, h_beliefs_p);
    hold on;
    params.schedule = 'sequential';
    figure(h_beliefs_s);
    plot_beliefs(H, q, num_points, params, h_beliefs_s);
    hold on;
end
figure(h_beliefs_p);
export_fig 'for_report/beliefs_p' '-pdf'
figure(h_beliefs_s)
export_fig 'for_report/beliefs_s' '-pdf'

num_points = 20;
% Шеннон в зависимости от r
n = 128;
q = 0.1;
ch = 1 + q * log2(q) + (1 - q) * log2(1 - q);
params.damping = 0.85;
k = [12:16:128];
j = 3;
err_bit = zeros(length(k), 1);
err_block = zeros(length(k), 1);
diver = zeros(length(k), 1);
for i = 1:length(k)
    display(['i --- ', num2str(i), ' / ', num2str(length(k))])
    H = make_ldpc_mex(n - k(i), n, j);
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

h = figure;
set(h, 'Color', 'w')
plot(k / n, err_bit, 'g')
hold on
plot(k / n, err_block, 'b')
plot(k / n, diver, 'r')
plot([ch, ch], [0, 1], 'm')
xlabel('r')
legend('bit error', 'block error', 'diver');
export_fig 'for_report/sh_r' '-pdf'
hold off

% Шеннон в зависимости от n
r = 0.25;
n = [16:16:128];
err_bit = zeros(length(n), 1);
err_block = zeros(length(n), 1);
diver = zeros(length(n), 1);
for i = 1:length(n)
    display(['i --- ', num2str(i), ' / ', num2str(length(n))])
    k = r * n(i);
    H = make_ldpc_mex(n(i) - k, n(i), j);
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

h = figure;
set(h, 'Color', 'w')
plot(n, err_bit, 'g')
hold on
plot(n, err_block, 'b')
plot(n, diver, 'r')
xlabel('n')
export_fig 'for_report/sh_n' '-pdf'
legend('bit error', 'block error', 'diver');
hold off

% Шеннон в зависимости от j
n = 128;
k = 32;
m = n - k;
j = [3:1:12];
err_bit = zeros(length(j), 1);
err_block = zeros(length(j), 1);
diver = zeros(length(j), 1);
for i = 1:length(j)
    display(['i --- ', num2str(i), ' / ', num2str(length(j))])
    H = make_ldpc_mex(m, n, j(i));
    [err_bit(i), err_block(i), diver(i)] = ldpc_mc(H, q, num_points, params);
end

h = figure;
set(h, 'Color', 'w')
plot(j, err_bit, 'g')
hold on
plot(j, err_block, 'b')
plot(j, diver, 'r')
xlabel('j')
legend('bit error', 'block error', 'diver');
export_fig 'for_report/sh_j' '-pdf'
hold off
