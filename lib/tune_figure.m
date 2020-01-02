% Fine tune the appearance of the figures

% Configuration
BoxLineWidth = 2;
AspectRatio = 1.26;
XLableDistance_y = - 0.09;     % Units = 'normalized'
YLableDistance_x = - 0.10;      % Units = 'normalized'
CanvasWidth = 560;              % Units = 'pixel'
CanvasHeight = 430;             % Units = 'pixel'
LineWidth = 2;
MarkerSize = 9;
MyFontSize = 14;

% Box
box on
set(gca, 'LineWidth', BoxLineWidth);

% Aspect ratio
pbaspect([AspectRatio, 1, 1]);

% Distance of labels
xlab = get(gca, 'XLabel');
set(xlab, 'Units', 'normalized');
set(xlab, 'Position', [0.5 XLableDistance_y 0]);

ylab = get(gca, 'YLabel');
set(ylab, 'Units', 'normalized');
set(ylab, 'Position', [YLableDistance_x 0.5 0]);

% Legend box
lgn = legend;
set(lgn, 'LineWidth', 1.5);

% Canvas size
set(gcf,'units','pixel');
gcf_pos = get(gcf, 'Position');
set(gcf,'Position', [gcf_pos(1:2) CanvasWidth CanvasHeight]);

% Change lines
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', LineWidth)
set(hline, 'MarkerSize', MarkerSize);

% Font size
set(gca, 'FontSize', MyFontSize);