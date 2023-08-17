% Script to plot previously generated data and reproduce Fig.2, 3 of
% "Foresight and relaxation enable efficient control of nonlinear complex systems", 
% by Xin Zhang, Xiaozhu Zhang, Gang Yan and Jack Murdoch Moore,
% currently under review at Physical Review Reesearch
%
% Jack Moore, August 2023
% 


close all;

% First produce Fig. 3:

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

ax.Units = 'Normalized';

axisInset = make_inset(patienceFig, ax, zoomArea, insetPos, insetLineCell);

exportgraphics(patienceFig, ['comparison-', vector_to_string(rList), '_', sysStr, '_cont-', vector_to_string(controlNodes), '_strat-', vector_to_string(strategyList), '.png']);


% Now produce Fig. 2:

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

%% Functions
%%%%%%%%%%%%%%%%%
%
% A function to turn a vector into a string format useful for file names
% 
% Jack Murdoch Moore, June 2022
% 
function str = vector_to_string(vec)
str = num2str(vec(1));
for ii = 2:numel(vec)
    str = [str, ',', num2str(vec(ii))];
end
end

function iaxes = make_inset(h, ax, ZoomArea, InsetPos, Lines)
    % Slightly modified (by Jack Murdoch Moore, 2023.02.19) from the following:
    % 
    % ************ MagInset.m *************
    % * Created: January 1, 2015          *
    % * By: Damien Frost                  *
    % * Inspiration By: Robert Richardson *
    % *************************************
    % 
    % MagInset (Magnified Inset) places an axes inset into an existing 
    % figure with handle h. Re-size your figure and axes before running
    % MagInset.
    %
    % h - handle to a figure on which to place an inset
    %
    % ax - handle to an axes object. Set to -1 to use the default axes in the
    %      figure.
    %
    % ZoomArea - a 4x1 array which defines the zoom area that inset
    %            should have. It is of the form: 
    %                [xmin xmax ymin ymax]
    %            The values are NOT normalized values, but graph values. In
    %            this way, you can easily read off the zoom area you would like
    %            to specify.
    %
    % InsetPos - defines the position of the inset. It is of the form:
    %                [xmin xmax ymin ymax]
    %            The values are NOT normalized values, but graph values. In
    %            this way, you can easily read off the inset area you would
    %            like to specify
    %
    % Lines - defines a list of lines to connect the zoom area with the inset
    %         graph. It can be empty. It should be of the form:
    %             Lines = {'NW','SW'; 'NE','SE'};
    %         It can have as many rows as you wish
    %         The first colum is the corner of the Zoom Area. The second column is the
    %         corner of the Inset Axes
    insetCol = 0.5*[1, 1, 1];
    lineWidth = 1;
    BadInput = 0;
    axesObjs = get(h, 'Children');  %axes handles
    % Determine which axes the inset will be placed on:
    if(ax == -1)
        MainAx = axesObjs(end);
    else
        MainAx = -1;
        for ii=1:1:max(size(axesObjs))
            if(axesObjs(ii) == ax)
                MainAx = axesObjs(ii);
                break;
            end
        end
        if(MainAx == -1)
            % Could not find the desired axes:
            fprintf('\nMagInset Error: Could not find the desired axes in the figure.\n');
            BadInput = 1;
        end
    end
    if(BadInput == 0)
        % Get the plot data:
        dataObjs = get(MainAx, 'Children');
        % Annotation positions are of the form:
        % [x y length height]
        % And are normalized to the figure
        % Calculate the normalize rectangular coordinates for the zoom area:
        [zax, zay] = xy2norm(MainAx, ZoomArea(1:2), ZoomArea(3:4));
        % Create the rectangle around the area we are going to zoom into:
        annotation('rectangle',[zax(1) zay(1) (zax(2) - zax(1)) (zay(2) - zay(1))], 'LineWidth', lineWidth, 'EdgeColor', insetCol);
        % Calculate the inset position in normalized coordinates;
        [ipx, ipy] = xy2norm(MainAx, InsetPos(1:2), InsetPos(3:4));
        if(nargin > 4)
            % Add the lines from the zoom area to the inset:
            numLine = size(Lines,1);
            if((numLine>0) && (size(Lines,2) == 2))
                lx = zeros(2,1);
                ly = zeros(2,1);
                for ii=1:1:numLine
                    jj = 1;
                    % Find the co-ordinate in the zoom area:
                    % y co-ordinates:
                    if(Lines{ii,jj}(1) == 'S')
                        ly(jj) = zay(1);
                    else
                        ly(jj) = zay(2);
                    end
                    % x co-ordinates:
                    if(Lines{ii,jj}(2) == 'W')
                        lx(jj) = zax(1);
                    else
                        lx(jj) = zax(2);
                    end
                    jj = 2;
                    % Find the co-ordinate in the inset axes:
                    % y co-ordinates:
                    if(Lines{ii,jj}(1) == 'S')
                        ly(jj) = ipy(1);
                    else
                        ly(jj) = ipy(2);
                    end
                    % x co-ordinates:
                    if(Lines{ii,jj}(2) == 'W')
                        lx(jj) = ipx(1);
                    else
                        lx(jj) = ipx(2);
                    end
                    % Add the line:
                    annotation('line', lx, ly, 'LineWidth', lineWidth, 'Color', insetCol);
                end
            end
        end
        % Add the second set of axes on the same plot:
        iaxes = axes('position', [ipx(1) ipy(1) (ipx(2) - ipx(1)) (ipy(2) - ipy(1))], 'LineWidth', lineWidth);
        hold on;
        box on;
        % Add the plots from the original axes onto the inset axes:
        copyobj(dataObjs,iaxes);
        iaxes.XTick = []; iaxes.YTick = [];
        iaxes.XColor = insetCol; iaxes.YColor = insetCol;
        % set the limits on the new axes:
        xlim(ZoomArea(1:2));
        ylim(ZoomArea(3:4));
        % Our work here is done.
    end
end
function [xn, yn] = xy2norm(axh, x, y)
    % ********* xy2norm.m *********
    % * Created: January 1, 2015  *
    % * By: Damien Frost          *
    % *****************************
    % This function takes a point (x, y) and returns the normalized
    % co-ordinates (xn, yn) given the axes with handle axh.
    %
    % *** Modifications: ***
    % 1) Added 'axh = handle(axh);' line to make the script backwards
    % compatible with MATLAB R2013a (Perhaps further back) - Thanks Shi
    % Zhao!
    axh = handle(axh);
    % Save some temporary variables about the axes:
    axPos = axh.Position;
    xlims = axh.XLim;
    ylims = axh.YLim;
    % Calculate the normalized co-ordinates:
    xn = axPos(1) + axPos(3) .* (x-xlims(1))./(xlims(2) - xlims(1));
    yn = axPos(2) + axPos(4) .* (y-ylims(1))./(ylims(2) - ylims(1));
    % GTFO:
end

%%%%%%%%%%%%%%%%%