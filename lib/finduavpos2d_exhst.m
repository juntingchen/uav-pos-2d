function [Fmin, Xopt] = finduavpos2d_exhst(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap)

BldLines = urbanMap.BldLines;
BldHeight = urbanMap.BldHeight;
BldTypes = urbanMap.BldTypes;

L = norm(PosBS(1:2) - PosUE(1:2), 2);
midpos = (PosBS(1:2) + PosUE(1:2)) / 2;
Xbd = midpos(1) + [-L/2, L/2];
Ybd = midpos(2) + [-L/2, L/2];
Xrange = Xbd(1): stepSizeMeter: Xbd(2);
Yrange = Ybd(1): stepSizeMeter: Ybd(2);

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
% h = U.Hdrone;
Hdrone = U.Hmin; % <- Note: we enforce the search plane to be at the minimum 
h = Hdrone;      %    height, just to reduce the search space in 2D.
            
Hb = U.Hbs;

Nx = length(Xrange);
Ny = length(Yrange);

Fmin = inf;
Xhat = [0, 0, 0];
for i = 1:Nx
    for j = 1:Ny
        Xpos = [Xrange(i), Yrange(j)];
        rhob = norm(Xpos - PosBS(1:2), 2);
        if rhob > L
            continue
        end
        
        rho = norm(Xpos(1:2) - PosUE(1:2), 2);
        Posdelta = Xpos(1:2) - PosUE(1:2);
        gamma= Posdelta(:).' * Uvec(:) / rho;
        theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));
        Cond = rho - L * cos(theta);
        if Cond > 0
            continue
        end
        
        los = IsLosK(PosUE, Xpos, BldLines, BldHeight, Hdrone, BldTypes);
        f = getcostf2DK([Xpos, Hdrone], [PosUE, 0], [PosBS, U.Hbs], los, U, fun);
        
        if f < Fmin
            Fmin = f;
            Xhat = Xpos;
        end
    end
end

Xopt = Xhat;