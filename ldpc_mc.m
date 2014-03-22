% Оценка характеристик LDPC-кода с помощью метода Монте Карло
%
% Вход:
%  H — проверочная матрица чётности, бинарная матрица размера m x n;
%  q — вероятность инверсии бита при передаче по каналу связи, число от 0 до 0.5;
%  num_points — общее количество экспериментов, число;
%
% Выход:
%  err_bit — вероятность битовой ошибки декодирования (относительно n бит 
%   кодового слова), число от 0 до 1, вычисляется по тем ситуациям, когда
%   алгоритм декодирования сошёлся (status < 2);
%  err_block — вероятность блоковой ошибки декодирования, число от 0 до 1, 
%   вычисляется по тем ситуациям, когда алгоритм декодирования сошёлся 
%   (status < 2)
%  diver — доля ситуаций расходимости алгоритма декодирования, число от 0 до 1.

function [err_bit, err_block, diver] = ldpc_mc(H, q, num_points)
    [m, n] = size(H);
    k = n - m;
    %display(['Channel ', num2str(1 + q * log2(q) + (1 - q) * log2(1 - q))]);
    %display(['Speed ', num2str(k / n)]);
    err_bit_num = zeros(n, 1);
    err_block_num = 0;
    diver_num = 0;
    good_status_num = 0;
    bad_status_num = 0;
    for iter = 1:num_points
        display(['#', num2str(iter)]);
        e = mod(binornd([1:n]', q), 2);
        s = mod(H * e, 2);
        [e_n, status] = ldpc_decoding(s, H, q, 'schedule', 'parallel', 'display', false);
        if status < 2
            good_status_num = good_status_num + 1;
            err_bit_num = err_bit_num + xor(e, e_n);
            err_block_num = err_block_num + any(e ~= e_n);
        else
            bad_status_num = bad_status_num + 1;
        end
    end
    err_bit = mean(err_bit_num) / good_status_num;
    err_block = err_block_num / good_status_num;
    diver = bad_status_num / num_points;
end
