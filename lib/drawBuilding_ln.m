function drawBuilding_ln(Blds, U)
% Draw buildings with single color filled with patterns. 
% Only for convex building shapes

if nargin < 2
    st = 40;
    LengthX = 1000;
    LengthY = 1000;
else
    st = U.BldSizeMin / 2;
    LengthX = U.LengthX;
    LengthY = U.LengthY;
    
end

Nbld = size(Blds, 1);
for ib = 1:Nbld
    Nlines = size(Blds{ib, 1}, 1);
    Lines = zeros(Nlines, 4);
    for j = 1:Nlines - 1
         Lines(j, :) = [Blds{ib, 1}(j, 1), Blds{ib, 1}(j, 2), ...
                        Blds{ib, 1}(j + 1, 1), Blds{ib, 1}(j + 1, 2)];
         % plot(Lines(j, [1 3]), Lines(j, [2 4]), 'b-', 'LineWidth', 2);
         H = line(Lines(j, [1 3]), Lines(j, [2 4]));
         set(H, 'color', 'b');
         set(H, 'LineStyle', '-');
         set(H, 'LineWidth', 2);
         if j == 1
             hold on
         end
    end
    Lines(Nlines, :) = [Blds{ib, 1}(Nlines, 1), Blds{ib, 1}(Nlines, 2), ...
                        Blds{ib, 1}(1, 1), Blds{ib, 1}(1, 2)];
    H = line(Lines(Nlines, [1 3]), Lines(Nlines, [2 4]));
    set(H, 'color', 'b');
    set(H, 'LineStyle', '-');
    set(H, 'LineWidth', 2);
    hold off
    
    % Fill pattern
    L = min(Blds{ib, 1}(:, 1));
    R = max(Blds{ib, 1}(:, 1));
    B = min(Blds{ib, 1}(:, 2));
    T = max(Blds{ib, 1}(:, 2));
    
    bH = T - B;
    bW = R - L;
    D = max(bH, bW);
    
    % st = 40;
    
    for x1 = L:st:(R + bH)
        y1 = B;
        x2 = x1 - D;
        y2 = y1 + D;
        XY1 = [x1, y1, x2, y2];
        
        T = lineSegmentIntersect(XY1, Lines);
        I = find(T.intAdjacencyMatrix(1, :) > 0);
        
        hold on
        if length(I) > 1
            H = line(T.intMatrixX(1, I(1:2)), T.intMatrixY(1, I(1:2)));
            set(H, 'color', 'b');
            set(H, 'LineStyle', '--');
            set(H, 'LineWidth', 1);
        end
        hold off
    end
    
    hold on
    xlim([0, LengthX]);
    ylim([0, LengthY]);
end

hold off