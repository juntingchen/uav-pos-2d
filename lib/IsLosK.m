function t = IsLosK(p1, p2, BldLines, BldHeight, DroneHeight, BldType)
% Version K: For two propagation segments. It returns values = 1 -
% ObstacleTypeCoefficient: 0 - LOS, 1 - k/(K-1) for the (k+1)th propagation
% segment, where virtual obstacle types = 0, 1, 2, ..., K - 1, and virtual
% obstacle coefficient = 0, 1/(K-1), 2/(K-1), ..., 1.
%
% NOT YET FINISHED!
%
% To test if there is a building that blocks the LOS between two points p1
% and p2 in 3D space.
%

Nbld = size(BldHeight, 1);

% We first project every thing on the XY plane to see whether the line
% segement intersect with any one of the buildings.
dU2D = norm(p1(1:2) - p2(1:2), 2);


ISLOS = 1;
for ib = 1:Nbld
    if inpoly(p1(1:2),  BldLines{ib}(:, 1:2))
        ISLOS = 0;
    else
        T = lineSegmentIntersect([p1(1:2), p2(1:2)], BldLines{ib});
        I = find(T.intAdjacencyMatrix(1, :) > 0);
        if ~isempty(I)
            BuildHeight = BldHeight(ib);
            for i_isec = 1:length(I)
                xI = T.intMatrixX(1, I(i_isec));    % The i_isec (th) intersec point
                yI = T.intMatrixY(1, I(i_isec));    % The i_isec (th) intersec point
                dU2I = norm(p1(1:2) - [xI yI], 2);      % distance, user to intersec point
                iHeight = DroneHeight * dU2I / dU2D;    % intersec point height at LOS
                if iHeight < BuildHeight
                    ISLOS = min(ISLOS, 1 - BldType(ib));
                    if ISLOS < 1e-9
                        break
                    end
                end
            end
        end
    end
end

t = ISLOS;
