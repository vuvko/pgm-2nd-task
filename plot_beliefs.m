% функция строит графики стабилизировывшихся beliefs

function plot_beliefs(H, q, num_points, params, h)

if nargin < 4 || ~isfield(params, 'color')
    params.color = 'b';
end
if ~isfield(params, 'style')
    params.style = '-';
end
if nargin < 5
    h = figure;
    set(h, 'Color', 'w');
end

figure(h);
[~, ~, ~, b_stat] = ldpc_mc(H, q, num_points, params);
plot(b_stat, params.style, 'Color', params.color, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Part of beliefs stabilazed');

end