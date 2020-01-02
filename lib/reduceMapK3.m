function [Blds, cAngles, BldTypes] = reduceMapK3(PosUE, fullBlds, fullBldTypes, U)
% Version K3: K=3 propagation segment case, identifying two types of
% obstacles. The new input variable fullBldTypes is a vector taking values
% 1 or 2, indicating the two types of buidlings (1 - wood, 2 - concrete
% wall)
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
BldTypes = zeros(Nfbld, 1);
CriticalBldInd = zeros(Nfbld, 1);
cnt = 0;

L = (U.LengthX + U.LengthY) / 2 / 2;
StepLength = 10;
Nsteps = max(1, round(2 * pi * L / StepLength));
thetas = 2 * pi / Nsteps: 2 * pi / Nsteps : 2 * pi;

cAngles = zeros(1, Nsteps);

k1BldIds = find(fullBldTypes == 1);
k2BldIds = find(fullBldTypes == 2);

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
        
        % [~, I] = max(BldAng);
        [~, I10] = max(BldAng(k1BldIds));
        [~, I20] = max(BldAng(k2BldIds));
        I1 = k1BldIds(I10);
        I2 = k2BldIds(I20);
        
        if ~CriticalBldInd(I2)
            
            cnt = cnt + 1;
            Blds{cnt, 1} = fullBlds{I2, 1};
            Blds{cnt, 2} = fullBlds{I2, 2};
            BldTypes(cnt) = 2;
            
            CriticalBldInd(I2) = 1;
        end
        
        if ~CriticalBldInd(I1) && BldAng(I1) > BldAng(I2) % the wood wall is in front of the concrete wall
            
            cnt = cnt + 1;
            Blds{cnt, 1} = fullBlds{I1, 1};
            Blds{cnt, 2} = fullBlds{I1, 2};
            BldTypes(cnt) = 1;
            
            CriticalBldInd(I1) = 1;
        end
        
        
    end
end

Blds = Blds(1:cnt, :);
BldTypes = BldTypes(1:cnt);
