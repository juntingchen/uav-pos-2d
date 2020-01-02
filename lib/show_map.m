function show_map(Xvec, Yvec, Zmat, BldArea, fid)

    if nargin < 5
        figure,
    else
        if ~isempty(fid)
            figure(fid),
        end
    end

    Mcolors = 64;
    cmap = jet(Mcolors);
    cmap(1, :) = [1 1 1];
    colormap(cmap(1:57,:));
    
    [Xmeter, Ymeter] = meshgrid(Xvec, Yvec);
    surf(Xmeter.', Ymeter.', Zmat, 'FaceColor','interp', 'edgecolor', 'none'); 
    if nargin >= 4 && ~isempty(BldArea)
        hold on
        for i = 1:size(BldArea, 1)
            bldx = BldArea{i, 1}; bldy = BldArea{i, 2};
            plot3(bldx, bldy, 50 * ones(size(bldy)), 'color', 'black');
        end
        hold off
    end
    % view(10, 5);
    view(0, 90);
    axis square
    xlim([Xvec(1), Xvec(end)]);
    ylim([Yvec(1), Yvec(end)]);
    set(gca, 'FontSize', 14);
    xlabel('X-axis [meter]');
    ylabel('Y-axis [meter]');

    % Tune the figure
    XLableDistance_y = - 0.08;     % Units = 'normalized'
    YLableDistance_x = - 0.10;      % Units = 'normalized'
    xlab = get(gca, 'XLabel');
    set(xlab, 'Units', 'normalized');
    set(xlab, 'Position', [0.5 XLableDistance_y 0]);
    ylab = get(gca, 'YLabel');
    set(ylab, 'Units', 'normalized');
    set(ylab, 'Position', [YLableDistance_x 0.5 0]);

%     % For the Washinton DC map
%     box on; BoxLineWidth = 2;
%     set(gca, 'LineWidth', BoxLineWidth);
%     xlim([0 800]);ylim([0 800]);
%     set(gca, 'XTick', 0:200:800);set(gca, 'YTick', 0:200:800);
%     
%     h = colorbar;
%     ylabel(h, 'Height [meter]', 'FontSize', 14);
%     set(h, 'YTick', 0:15:45); 
end