function [Blds, cAngles] = reduceMap(PosUE, fullBlds, U)
% Version 2: added the statsicis, i.e., cAngles (for 'criticalAngles'), Prob_LOS(elevation_angle)

DroneHeight = U.Hdrone;

Nfbld = size(fullBlds, 1);

BldLines = cell(Nfbld, 1);
for ib = 1:Nfbld
    Nlines = size(fullBlds{ib, 1}, 1);
    BldLines{ib} = zeros(Nlines, 4);
    for j = 1:Nlines - 1
        BldLines{ib}(j, :) = [fullBlds{ib, 1}(j, 1), fullBlds{ib, 1}(j, 2), ...
                              fullBlds{ib, 1}(j + 1, 1), fullBlds{ib, 1}(j + 1, 2)];
    end
    BldLines{ib}(Nlines, :) = [fullBlds{ib, 1}(Nlines, 1), fullBlds{ib, 1}(Nlines, 2), ...
                              fullBlds{ib, 1}(1, 1), fullBlds{ib, 1}(1, 2)];
end

Blds = cell(Nfbld, 2);
CriticalBldInd = zeros(Nfbld, 1);
cnt = 0;

L = (U.LengthX + U.LengthY) / 2 / 2;
StepLength = 10;
Nsteps = max(1, round(2 * pi * L / StepLength));
thetas = 2 * pi / Nsteps: 2 * pi / Nsteps : 2 * pi;

cAngles = zeros(1, Nsteps);

for i = 1:Nsteps
    ta = thetas(i);
    PosD = L * [cos(ta), sin(ta)] + PosUE(1:2);
    
    XYu = PosUE;
    XYd = PosD;

    dU2D = norm(XYu - XYd, 2);      % distance, user to drone (ground projected)
    BldAng = zeros(Nfbld, 1);
    
    ISLOS = 1;
    for ib = 1:Nfbld
        T = lineSegmentIntersect([XYu XYd], BldLines{ib});
        I = find(T.intAdjacencyMatrix(1, :) > 0);
        if ~isempty(I)
            BuildHeight = fullBlds{ib, 2};
            angle = 0;
            for i_isec = 1:length(I)
                xI = T.intMatrixX(1, I(i_isec));    % The i_isec (th) intersec point
                yI = T.intMatrixY(1, I(i_isec));    % The i_isec (th) intersec point
                dU2I = norm(XYu - [xI yI], 2);      % distance, user to intersec point
                iHeight = DroneHeight * dU2I / dU2D;    % intersec point height at LOS
                if iHeight < BuildHeight
                    % there is intersection
                    ISLOS = 0;
                    ang = atan(BuildHeight / dU2I);
                    if ang > angle
                        angle = ang;
                    end
                    % break;
                end
            end
            BldAng(ib) = angle;
        end
    end

    if ISLOS
        cAngles(i) = 0;
        
    elseif ~ISLOS
        cAngles(i) = angle;
        
        [~, I] = max(BldAng);
        
        if ~CriticalBldInd(I)
            
            cnt = cnt + 1;
            Blds{cnt, 1} = fullBlds{I, 1};
            Blds{cnt, 2} = fullBlds{I, 2};
            
            CriticalBldInd(I) = 1;
        end
    end
end

Blds = Blds(1:cnt, :);
