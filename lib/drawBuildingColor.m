function drawBuildingColor(Blds, U)
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

maxH = 0;
for ib = 1:Nbld
    H = Blds{ib, 2};
    if H > maxH
        maxH = H;
    end
end
Mcolors = ceil(maxH);
mycolors = jet(Mcolors);

for ib = 1:Nbld
    
    VX = Blds{ib, 1}(:, 1);
    VY = Blds{ib, 1}(:, 2);
    H = Blds{ib, 2};
    fill(VX, VY, mycolors(min(Mcolors, max(1, round(H))), :));
    
%     Nlines = size(Blds{ib, 1}, 1);
%     Lines = zeros(Nlines, 4);
%     for j = 1:Nlines - 1
%          Lines(j, :) = [Blds{ib, 1}(j, 1), Blds{ib, 1}(j, 2), ...
%                         Blds{ib, 1}(j + 1, 1), Blds{ib, 1}(j + 1, 2)];
%          plot(Lines(j, [1 3]), Lines(j, [2 4]), 'b-', 'LineWidth', 2);
%          if j == 1
%              hold on
%          end
%     end
%     Lines(Nlines, :) = [Blds{ib, 1}(Nlines, 1), Blds{ib, 1}(Nlines, 2), ...
%                         Blds{ib, 1}(1, 1), Blds{ib, 1}(1, 2)];
%     plot(Lines(Nlines, [1 3]), Lines(Nlines, [2 4]), 'b-', 'LineWidth', 2);

    
    hold on
    xlim([0, LengthX]);
    ylim([0, LengthY]);
end

hold off

% colormap(mycolors);

