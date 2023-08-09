% Script to plot previously generated data and reproduce Fig. 4 of
% "Foresight and relaxation enable efficient control of nonlinear complex systems", 
% by Xin Zhang, Xiaozhu Zhang, Gang Yan and Jack Murdoch Moore,
% currently under review at Physical Review Reesearch
%
% Jack Moore, August 2023
% 

close all;
clear;

load('comparison-1e-06,0.001,0.1,1_bist_cont-3,5_strat-43,43,43,43_energy-1,1,1,1.mat');

plotDim1 = 1;
plotDim2 = 3;
plotDim3 = 5;

num_r = numel(rList);

patienceFig = figure; hold on;
patienceLegCell = {};

colOrder = get(gca, 'ColorOrder');
colOrder([3, 4], :) = colOrder([4, 3], :);
colOrder = [colOrder; colOrder];
lineStyleCell = {'-', '-', '--', '-.', ':'};
lineStyleCell = [lineStyleCell, lineStyleCell];
lineWidth = 1;
for ii_r = 1:num_r
    plot(NaN, NaN, 'LineWidth', 2*lineWidth, 'Color', 0*colOrder(ii_r + 1, :), 'LineStyle', lineStyleCell{ii_r + 1});
end

fontSize = 15;
axis square;
axis tight;
box on;
set(gca, 'FontSize', fontSize, 'TickLabelInterpreter', 'LaTeX');
plot([x0(plotDim1), xf(plotDim1)], [x0(plotDim2), xf(plotDim2)], 'ko', 'LineWidth', 2*lineWidth);

for ii_r = 1:num_r
    r = rList(ii_r);
    
    rStr = num2str(r);
    if contains(rStr, 'e')
        rStr = [strrep(rStr, '1e', '10^{'), '}'];
    end
    
    patienceLegCell = [patienceLegCell, ['$r=', rStr, '$']];

    x_total = x_totalCell{ii_r};
    t_total = t_totalCell{ii_r};
    d_total = d_totalCell{ii_r};
    s_total = s_totalCell{ii_r};
    e_total = e_totalCell{ii_r};
    
    x = x_total(plotDim1, :);
    y = x_total(plotDim2, :);
    z = x_total(plotDim3, :);
    % c = log10(e_total);
    % c = [diff(e_total)./diff(t_total), NaN];
    c = log10(e_total); propertyNameStr = '$\log_{10}$Energy';
    % c = e_total; propertyNameStr = 'Energy';
    % c = e_total/max(e_total); propertyNameStr = 'Normalised energy';
    c(c < 1) = 1;
    
surface([x; x], [y; y], [z; z], [c; c],...
        'FaceCol','None',...
        'EdgeCol','Interp',...
        'LineWidth', 2*lineWidth,...
        'LineStyle', lineStyleCell{ii_r + 1});
    
end

cb = colorbar; colormap(jet);
cb.TickLabelInterpreter = 'LaTeX';
cb.Label.Interpreter = 'LaTeX';
cb.Label.String = propertyNameStr;

xlabel(['$x_', num2str(plotDim1), '$'], 'Interpreter', 'LaTeX'); ylabel(['$x_', num2str(plotDim2), '$'], 'Interpreter', 'LaTeX'); zlabel(['$x_', num2str(plotDim3), '$'], 'Interpreter', 'LaTeX');

text(x0(plotDim1), x0(plotDim2), '${\textbf{\textit{x}}}_{{I}}$', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom', 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Margin', eps);
text(xf(plotDim1), xf(plotDim2), '${\textbf{\textit{x}}}_{{F}}$', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top', 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Margin', eps);
leg = legend(patienceLegCell, 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'Box', 'Off');
legPos = leg.Position;
legPos(2) = 0.15;
leg.Position = legPos;

exportgraphics(patienceFig, ['comparison-', vector_to_string(rList), '_', sysStr, '_cont-', vector_to_string(controlNodes), '_strat-', vector_to_string(strategyList), '.png']);