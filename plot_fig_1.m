load './data/fig_1_2x24x4_x.mat'
load './data/fig_1_2x24x4_y.mat'

h = figure();
semilogy(xdata{4}, ydata{4}, '-b *', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'MarkerSize', 7);
hold on;
semilogy(xdata{2}, ydata{2}, '-g o', 'LineWidth', 1.5, 'Color', "#0072BD", 'MarkerFaceColor', "#FFFF00");
semilogy(xdata{1}, ydata{1}, '-r d', 'LineWidth', 1.5, 'Color', "#FF0000", 'MarkerFaceColor', "#FFFF00");


load './data/fig_2_2x24x4_6_x.mat'
load './data/fig_2_2x24x4_6_y.mat'
semilogy(xdata{3}, ydata{3}, '--*', 'LineWidth', 1.5, 'Color', "#7E2F8E", 'MarkerSize', 7);
semilogy(xdata{2}, ydata{2}, '--g o', 'LineWidth', 1.5, 'Color', "#0072BD", 'MarkerFaceColor', "#FFFF00");
semilogy(xdata{1}, ydata{1}, '--r d', 'LineWidth', 1.5, 'Color', "#FF0000", 'MarkerFaceColor', "#FFFF00");


grid minor;
ylabel('CRB (dB)', 'FontSize', 14, 'Interpreter','latex');
xlabel('SNR (dB)', 'FontSize', 14, 'Interpreter','latex');
legend('Unstructured: OP', 'Structured ULA: OP', 'Structured UCyA: OP', ...
    'Unstructured: SB', 'Structured ULA: SB', 'Structured UCyA: SB', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Edgecolor', 'white');
hAx=gca;                              % get the axes handle
hAx.XTickLabel=hAx.XTickLabel;        % overwrite the existing tick labels with present values
set(gcf,'color','w');
xticks('auto');
set(gca,'FontName','Times','fontsize', 12);