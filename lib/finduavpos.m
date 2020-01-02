function [Fmin, Xopt, stepCount, XposArray2d] = finduavpos(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap)

BldLines = urbanMap.BldLines;
BldHeight = urbanMap.BldHeight;
BldTypes = urbanMap.BldTypes;

PosUE3 = [PosUE, 0];
PosBS3 = [PosBS, U.Hbs];

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
maxStep = ceil((2.4 * U.K - 1.4) * norm(PosBS(1:2) - PosUE(1:2), 2) / stepSizeMeter);

Fmin = inf;
Xopt = [0, 0];

XposArray2d = zeros(maxStep, 2);
FArray2d = zeros(maxStep, 1);

Xpos2d = PosBS;
Stage = 1;
cnt2d = 0;
Rhos2d = zeros(1, maxStep); % Record for Stage 1, rho value along the way
ks = zeros(1, maxStep);   % Record for stage 1, segment label along the way
ksearch = 1;
while Stage < 4 && cnt2d < maxStep
    cnt2d = cnt2d + 1;
    
    Xpos0 = Xpos2d;   % UAV Position
    
    los = IsLosK(PosUE, Xpos0, BldLines, BldHeight, U.Hdrone, BldTypes);
    ksegment = round((1 - los) * (U.K - 1) + 1);   % propagation segment index, k = 1,2,...,K
    f = getcostf2DK([Xpos0, U.Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
    FArray2d(cnt2d) = f;
    XposArray2d(cnt2d, :) = Xpos0;
    if f < Fmin
        Fmin = f;
        Xopt = Xpos0;
    end
    
    if Stage == 1
        % Stage 1: Search on the User-BS axis
        if ksegment == 1 % LOS segment
            Rhos2d(cnt2d) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt2d) = 1;
            Rhos02d = findCriticalPoints2DK(Rhos2d, ks, [Xpos0, U.Hdrone], PosUE3, PosBS3, U, fun);
            
            rho0 = Rhos02d(ksearch);
            theta0 = stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
            
            Stage = 2;
            
        elseif norm(PosUE - Xpos0) > stepSizeMeter 
            Rhos2d(cnt2d) = norm(PosUE(1:2) - Xpos0(1:2), 2);
            ks(cnt2d) = ksegment;
            searchDirection = (PosUE - PosBS) / norm(PosUE - PosBS);
            Xpos2d = Xpos0 + searchDirection * stepSizeMeter;
            
        else
            % Indoor case
            Stage = 4;
        end
        
    elseif Stage == 2
        % Stage 2: Search on the right branch
        searchDirection = uavSearchDirection2DK([Xpos2d, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos2d = Xpos0 + searchDirection * stepSizeMeter;
        
        if ksegment <= ksearch % LOS or virtual LOS
            uavColor = [1, 0, 0];
        else
            uavColor = [0, 1, 0];
        end
        
        if norm(searchDirection) < 1e-10
            rho0 = Rhos02d(ksearch);
            theta0 = - stepSizeMeter / rho0;
            M = [cos(theta0), -sin(theta0)
                 sin(theta0), cos(theta0)];
            Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
            
            Stage = 3;
        end
        
    elseif Stage == 3
        % Stage 3: Search on the left branch
        searchDirection = uavSearchDirection2DK([Xpos2d, U.Hdrone], ...
                                               PosUE3, PosBS3, ksegment, ksearch, U, fun);
        %
        Xpos2d = Xpos0 + searchDirection * stepSizeMeter;

        if norm(searchDirection) < 1e-10
            
            ksearch = ksearch + 1;
            if ksearch < U.K
                rho0 = Rhos02d(ksearch);
                theta0 = stepSizeMeter / rho0;
                M = [cos(theta0), -sin(theta0)
                     sin(theta0), cos(theta0)];
                Xpos2d = PosUE + (rho0 * M * Uvec(:)).';
                Stage = 2;
            else
                Stage = 4; % Algorithm terminates
            end
        end
        
    else
        % The entire algorithm terminates
    end
end

stepCount = cnt2d;
XposArray2d = XposArray2d(1:stepCount, :);