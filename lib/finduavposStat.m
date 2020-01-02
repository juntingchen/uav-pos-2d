function [Fmin, Xopt, stepCount] = finduavposStat(PosUE, PosBS, U, fun, stepSizeMeter, urbanMap, losStat)
% Statistical optimization baseline

BldLines = urbanMap.BldLines;
BldHeight = urbanMap.BldHeight;
% BldTypes = urbanMap.BldTypes;
% 
% PosUE3 = [PosUE, 0];
% PosBS3 = [PosBS, U.Hbs];

Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);

% Exhaustive algorithm (to be optimized)
L = norm(PosBS(1:2) - PosUE(1:2), 2);
rho_array = L - abs(stepSizeMeter) : - abs(stepSizeMeter): abs(stepSizeMeter);
Nsteps = length(rho_array);
Fvec = zeros(1, Nsteps);
for i = 1:Nsteps
    rho = rho_array(i);
    elev_angle = atan(U.Hdrone / rho);
    
    [~, I] = min(abs(losStat.ElvAngles - elev_angle));
    P_Los = losStat.Plos(I); 
    P_NLos = 1 - P_Los;
    
    d1 = sqrt(rho ^ 2 + (U.Hdrone - U.Huser) ^ 2);
    gDrone = P_Los * 10 .^ ((log10(d1) * U.A1 + U.B1) / 10) / U.Noise ...
            + P_NLos * 10 .^ ((log10(d1) * U.A2 + U.B2) / 10) / U.Noise;
        
    d0 = sqrt((L - rho)^2 + (U.Hdrone - U.Hbs)^2);
    gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
    
    Fvec(i) = fun(gDrone, gBS); 
end
[~, I] = min(Fvec);
rho_opt = rho_array(I);

Xopt = PosUE + (rho_opt * Uvec(:)).';

% Evaluation on actual environment
if IsLos(PosUE, Xopt, BldLines, BldHeight, U.Hdrone)
    % LOS case
    d1 = sqrt(rho_opt ^ 2 + (U.Hdrone - U.Huser) ^ 2);
    gDrone = 1 * 10 .^ ((log10(d1) * U.A1 + U.B1) / 10) / U.Noise ...
            + 0 * 10 .^ ((log10(d1) * U.A2 + U.B2) / 10) / U.Noise;
        
    d0 = sqrt((L - rho_opt)^2 + (U.Hdrone - U.Hbs)^2);
    gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
    
    Fmin = fun(gDrone, gBS); 
else
    % NLOS case
    d1 = sqrt(rho_opt ^ 2 + (U.Hdrone - U.Huser) ^ 2);
    gDrone = 0 * 10 .^ ((log10(d1) * U.A1 + U.B1) / 10) / U.Noise ...
            + 1 * 10 .^ ((log10(d1) * U.A2 + U.B2) / 10) / U.Noise;
        
    d0 = sqrt((L - rho_opt)^2 + (U.Hdrone - U.Hbs)^2);
    gBS = 10 .^ ((log10(d0) * U.A0 + U.B0) / 10) / U.Noise;
    
    Fmin = fun(gDrone, gBS); 
end

stepCount = Nsteps;
