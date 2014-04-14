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
