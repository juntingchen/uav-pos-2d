function f = getcostf2DK(PosUAV, PosUE, PosBS, los, U, fun)
% Version 2DK, for 2D UAV position, where the height is not optimized
% Version K, for K segmenet case, 
% where parameter los = 1, ..., 1/(K-1), 0
%
% Compute the partial optimal cost at a UAV position in 2D. 
% 
% INPUT 
%   PosUAV, PosUE, PosBS, positions in 3D
%   U, configuration structure, containing channel parameters
%   fun, objective function in terms of the channel gains
%
% OUTPUT
%   fmin,   the minimum cost over the ray with phi evalvation angle
%           depending on the UAV-user position
%   lopt,   corresponding coordinate in l.
%   Xopt    Optimal UAV position
%
% Constraint: l >= sqrt(rho^2 + h^2)
% NOTE: Assume user height = 0.

    L = norm(PosUE(1:2) - PosBS(1:2), 2);
    rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
    Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
    Posdelta = PosUAV(1:2) - PosUE(1:2);

    gamma= Posdelta(:).' * Uvec(:) / rho;
    theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

    h = PosUAV(3);
    Hb = PosBS(3);

    lmin = sqrt(rho^2 + h^2);
    % In the 2D case, we do not optimize the height of the UAV! %
    lopt = lmin;
    f = fl(lopt, rho, theta, L, h, Hb, U, fun, los);

end

function f = fl(l, rho, theta, L, h, Hb, U, fun, los)
    k = round((1 - los) * (U.K - 1) + 1);   % propagation segment index
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    gainu = 10 ^ ((U.Alpha(k) * log10(du2u) + U.Beta(k)) / 10) / U.Noise;

    f = fun(gainu, gainb);
end