% Script to plot previously generated data and reproduce Fig. 3 of
% "Foresight and relaxation enable efficient control of nonlinear complex systems", 
% by Xin Zhang, Xiaozhu Zhang, Gang Yan and Jack Murdoch Moore,
% currently under review at Physical Review Reesearch
%
% Jack Moore, August 2023
% 

close all;
clear;

load('comparison-0.0001,0.0001,0.0001_bist_cont-3,4_strat-45,46,43_energy-100000,1,1.mat');

plotDim1 = 1;
plotDim2 = 3;
plotDim3 = 5;

num_r = numel(rList);

zoomArea = [-4, -3.85, -0.6, 3.6];
insetPos = [0, 4, -1.5, 2.5];
insetLineCell = {'SW', 'SW'; 'NE', 'NE'};
% insetLineCell = {'SW', 'SW'};

patienceFig = figure; hold on;
patienceLegCell = {};

ax = gca;
colOrder = get(ax, 'ColorOrder');
colOrder([3, 4], :) = colOrder([4, 3], :);
colOrder = [colOrder; colOrder];
lineStyleCell = {'-', '-', '--', '-.', ':'};
lineStyleCell = [lineStyleCell, lineStyleCell];
lineWidth = 1;
for ii_q = 1:num_r
    plot(NaN, NaN, 'LineWidth', 2*lineWidth, 'Color', 0*colOrder(ii_q + 1, :), 'LineStyle', lineStyleCell{ii_q + 1});
end

fontSize = 15;
axis square;
% daspect([1, 1, 1]);
%title(['mode=', num2str(mo), ', p=', num2str(q), ', dt=', num2str(dt), ',dx=', num2str(dx), ',de=', num2str(de), ',T_0=', num2str(T0), ', ', num2str(control_nodes), ': ', num2str(dist0), '\rightarrow ', num2str(dist)]);
box on;
set(gca, 'FontSize', fontSize, 'TickLabelInterpreter', 'LaTeX');
plot([x0(plotDim1), xf(plotDim1)], [x0(plotDim2), xf(plotDim2)], 'ko', 'LineWidth', 2*lineWidth);

for ii_q = 1:num_r
    strategy = strategyList(ii_q);
    
    switch strategy
        case 45
            modeStr = 'EILOCS';
        case 46
            modeStr = 'DDLOCS';
        case 43
            modeStr = 'ALITE';
    end
    
    patienceLegCell = [patienceLegCell, modeStr];

    x_total = x_totalCell{ii_q};
    t_total = t_totalCell{ii_q};
    d_total = d_totalCell{ii_q};
    s_total = s_totalCell{ii_q};
    e_total = e_totalCell{ii_q};
    
    x = x_total(plotDim1, :);
    y = x_total(plotDim2, :);
    z = x_total(plotDim3, :);
    % c = log10(e_total);
    % c = [diff(e_total)./diff(t_total), NaN];
    c = log10(e_total); propertyNameStr = '$\log_{10}$Energy';
    % c = e_total; propertyNameStr = 'Energy';
    % c = e_total/max(e_total); propertyNameStr = 'Normalised energy';
    c(c < 2) = 2;
    
surface([x; x], [y; y], [z; z], [c; c],...
        'FaceCol','None',...
        'EdgeCol','Interp',...
        'LineWidth', 2*lineWidth,...
        'LineStyle', lineStyleCell{ii_q + 1});
    
end

ylim([-4, 6]);

cb = colorbar; colormap(jet);
cb.TickLabelInterpreter = 'LaTeX';
cb.Label.Interpreter = 'LaTeX';
cb.Label.String = propertyNameStr;

xlabel(['$x_', num2str(plotDim1), '$'], 'Interpreter', 'LaTeX'); ylabel(['$x_', num2str(plotDim2), '$'], 'Interpreter', 'LaTeX'); zlabel(['$x_', num2str(plotDim3), '$'], 'Interpreter', 'LaTeX');

text(x0(plotDim1), x0(plotDim2), '${\textbf{\textit{x}}}_{{I}}$', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom', 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Margin', eps);
text(xf(plotDim1), xf(plotDim2), '${\textbf{\textit{x}}}_{{F}}$', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top', 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Margin', eps);
leg = legend(patienceLegCell, 'Location', 'South', 'Interpreter', 'LaTeX', 'Box', 'Off');

axisInset = make_inset(patienceFig, ax, zoomArea, insetPos, insetLineCell);

exportgraphics(patienceFig, ['comparison-', vector_to_string(rList), '_', sysStr, '_cont-', vector_to_string(controlNodes), '_strat-', vector_to_string(strategyList), '.png']);