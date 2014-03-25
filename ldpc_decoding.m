% Алгоритм декодирования LDPC-кода в синдромном представлении
%
% Вход:
%  s — наблюдаемый синдром, бинарный вектор-столбец длины m;
%  H — проверочная матрица чётности, бинарная матрица размера m x n;
%  q — вероятность инверсии бита при передаче по каналу связи,
%   число от 0 до 0.5;
%  (param_name, param_value) — набор необязательных параметров алгоритма, 
%   следующие имена и значения возможны:
%    'schedule' — расписание пересчёта сообщений, возможные значения 
%     'parallel' и 'sequential', по умолчанию = 'parallel';
%    'damping' — коэффициент дэмпфирования при пересчёте сообщений, 
%     число от 0 до 1, по умолчанию = 1;
%    'max_iter' — максимальное число итераций алгоритма декодирования, 
%     число, по умолчанию = 200;
%    'eps' — порог стабилизации для beliefs, число, по умолчанию = 1e-4;
%    'display' — режим отображения, true или false, если true, 
%     то отображается промежуточная информация на итерациях, например, 
%     номер итерации, текущее число ошибок декодирования, невязка для
%     сообщений и т.д.
%
% Выход:
%  e — восстановленный вектор ошибок, бинарный вектор-столбец длины n;
%  status — результат декодирования: 
%    равен 0, если найден вектор e, соответствующий входному синдрому s,
%    равен 1, если произошла стабилизация значений beliefs, 
%    равен 2, если произошел выход по максимальному числу итераций.

function [e, status] = ldpc_decoding(s, H, q, varargin)

[m, n] = size(H);
e = zeros(n, 1);
status = -1;

%-------------------
% Разбор аргументов
[err, schedule, damping, max_iter, eps, dout] = parse_arg(varargin{:});
if err > 0
    v = varargin{err+1};
    if isnumeric(v)
        v = num2str(v);
    end
    error('LDPC_Decoding:UknownArgumentValue', ...
        'Argument %s has unknown value %s.', ...
        varargin{err}, v);
elseif err < 0
    error('LDPC_Decoding:UknownArgument', ...
        'Unknown argument %s.', varargin{-err});
end

%-------------------
% Инициализация

M_to_1 = zeros(size(H));
M_to_0 = zeros(size(H));
M_from_1 = zeros(size(H));
M_from_0 = zeros(size(H));
%e = mod(binornd([1:n]', q), 2);
b_1 = ones(n, 1);
b_0 = ones(n, 1);
% стартуем из унарных потенциалов
M_to_1 = M_to_1 + q;
M_to_0 = 1 - M_to_1;

% определяем последовательность обработки вершин/факторов
%sq = [[1:m], -[1:n]];
sq = [-[1:n], [1:m]];

%-------------------
% Основной цикл

for iter = 1:max_iter
    if dout
        display(['Iteration ', num2str(iter), ':']);
    end
    %-------------------
    % Пересчитываем сообщения вершин/факторов
    for ind = 1:length(sq)
        if sq(ind) > 0
            % пересчитываем сообщения от фактора к вершинам
            j = sq(ind);
            %delta = repmat(M_to_0(j, :) - M_to_1(j, :), n, 1);
            %idx = repmat(H(j, :), n, 1);
            %idx(eye(n) == 1) = 0;
            %delta(idx == 0) = 1;
            %dr = prod(delta, 2)';
            %if s(j) == 1
            %    M_from_0(j, :) = (1 - dr) / 2;
            %    M_from_1(j, :) = (1 + dr) / 2;
            %else
            %    M_from_0(j, :) = (1 + dr) / 2;
            %    M_from_1(j, :) = (1 - dr) / 2;               
            %end
            %% нормируем сообщения
            %%nm = (M_from_0(j, :) + M_from_1(j, :));
            %%M_from_0(j, :) = M_from_0(j, :) ./ nm;
            %%M_from_1(j, :) = M_from_1(j, :) ./ nm;
            
            f = M_from_0;
            idx = H;
            idx(j, :) = 0;
            f(idx == 0) = 1;
            M_to_0(j, :) = (1 - q) * prod(f, 1);
            f = M_from_1;
            f(idx == 0) = 1;
            M_to_1(j, :) = q * prod(f, 1);
            
            nm = M_to_0(j, :) + M_to_1(j, :);
            M_to_0(j, :) = M_to_0(j, :) ./ nm;
            M_to_1(j, :) = M_to_1(j, :) ./ nm;            
        else
            % пересчитываем сообщения от вершины к факторам
            i = -sq(ind);
            %f = repmat(M_from_0(:, i), 1, m);
            %idx = repmat(H(:, i), 1, m);
            %idx(eye(m) == 1) = 0;
            %f(idx == 0) = 1;
            %M_to_0(:, i) = (1 - q) * prod(f, 1)'; 
            %f = repmat(M_from_1(:, i), 1, m);
            %f(idx == 0) = 1;
            %M_to_1(:, i) = q * prod(f, 1)';
            %% нормируем сообщения
            %nm = (M_to_0(:, i) + M_to_1(:, i));
            %M_to_0(:, i) = M_to_0(:, i) ./ nm;
            %M_to_1(:, i) = M_to_1(:, i) ./ nm;
            
            delta = M_to_0 - M_to_1;
            idx = H;
            idx(:, i) = 0;
            delta(idx == 0) = 1;
            dr = (1 + prod(delta, 2)) / 2;
            dr(s == 1) = 1 - dr(s == 1);
            M_from_0(:, i) = dr;
            M_from_1(:, i) = 1 - dr;
        end
        %M_to_0
        %M_from_0
        %pause;
    end
    %-------------------
    % Вычисляем beliefs
    b_old = [b_0, b_1];
    dq = M_to_1;
    dq(H == 0) = 1;
    b_1 = q .* prod(dq, 1)';
    df = M_to_0;
    df(H == 0) = 1;
    b_0 = (1 - q) .* prod(df, 1)';
    %-------------------
    % Оценка вектора ошибок
    [~, e] = max([b_0, b_1], [], 2);
    e = e - 1;
    %-------------------
    % Проверка критериев остановки
    db = max(max(abs(b_old - [b_0, b_1])));
    if dout
        display(['  Difference in beliefs = ', num2str(db), ';']);
        display(['  Bit error = ', num2str(sum(xor(mod(H * e, 2), s)) / n), ';']);
    end
    if all(mod(H * e, 2) == s)
        if dout
            disp 'Exiting with status 0.';
        end
        status = 0;
        return;
    elseif db < eps
        if dout
            disp 'Exiting with status 1.';
        end
        status = 1;
        return;
    end
    % меняем последовательность обработки вершин/факторов
    % при последовательном расписании
    if strcmp(schedule, 'sequential')
        sq = sq(randperm(length(sq)));
    end
end
status = 2;
if dout
    disp 'Exiting with status 2.';
end

end

function [err, schedule, damping, max_iter, eps, dout] = parse_arg(varargin)

% Установка параметров по-умолчанию
schedule = 'parallel';
damping = 1;
max_iter = 200;
eps = 1e-4;
dout = true; % для постоянного дебага в процессе написания реализации.
err = 0;

for i = 1:2:length(varargin)
    if strcmp(varargin{i}, 'schedule')
        if strcmp(varargin{i+1}, 'sequential')
            schedule = 'sequential';
        elseif strcmp(varargin{i+1},'parallel')
            schedule = 'parallel';
        else
            err = i;
            return;
        end
    elseif strcmp(varargin{i}, 'damping')
        if varargin{i+1} > 0 && varargin{i+1} <= 1
            damping = varargin{i+1};
        else
            err = i;
            return;
        end
    elseif strcmp(varargin{i}, 'max_iter')
        if varargin{i+1} > 0
            max_iter = varargin{i+1};
        else
            err = i;
            return;
        end
    elseif strcmp(varargin{i}, 'eps')
        eps = varargin{i+1};
    elseif strcmp(varargin{i}, 'display')
        if strcmp(varargin{i+1}, 'true') || varargin{i+1} == true
            dout = true;
        elseif strcmp(varargin{i+1}, 'false') || varargin{i+1} == false
            dout = false;
        else
            err = i;
            return;
        end
    else
        err = -i;
        return;
    end
end

end
