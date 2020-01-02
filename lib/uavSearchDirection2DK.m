function D = uavSearchDirection2DK(PosUAV, PosUE, PosBS, ksegment, ksearch, U, fun)
% Version K, K segment case
% INPUT 
%   PosUAV, PosUE, PosBS, positions in 3D
%   U, configuration structure, containing channel parameters
%   fun, objective function in terms of the channel gains
%   stepSize, the step size

L = norm(PosUE(1:2) - PosBS(1:2), 2);
rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
Posdelta = PosUAV(1:2) - PosUE(1:2);

gamma= Posdelta(:).' * Uvec(:) / rho;
theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

h = PosUAV(3);
Hb = PosBS(3);

% ---- Numerical Computation of the Partial Derivatives ---- %
delta = 1e-4;

F_plus = F(rho + delta, theta, L, h, Hb, U, fun, ksearch);
F_minus = F(rho - delta, theta, L, h, Hb, U, fun, ksearch);
dF_drho = (F_plus - F_minus) / (2 * delta);

F_plus = F(rho, theta + delta, L, h, Hb, U, fun, ksearch);
F_minus = F(rho, theta - delta, L, h, Hb, U, fun, ksearch);
dF_dtheta = (F_plus - F_minus) / (2 * delta);


if dF_drho >= 0 || rho > L * cos(theta) % Check stopping criteria
    D = zeros(1, 2);
else
    if ksegment <= ksearch % LOS or virtual LOS
        D = Posdelta / norm(Posdelta, 2);
    else
        deltaX = Rot(theta) * Uvec(:) + rho * dRot(theta) * Uvec(:) ...
                 * ( - dF_drho / dF_dtheta);
        D = deltaX.' / norm(deltaX);
    end
        
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = Rot(ta)
    T = [cos(ta)    -sin(ta)
         sin(ta)    cos(ta)];
end

function T = dRot(ta)
    T = [-sin(ta)    -cos(ta)
         cos(ta)    -sin(ta)];
end



function f = F(rho, theta, L, h, Hb, U, fun, ksearch)
    % compute the function F(\rho,\theta) in the paper

    lopt = sqrt(rho^2 + h^2);
    f = fl(lopt, rho, theta, L, h, Hb, U, fun, ksearch);
    
end

function f = fl(l, rho, theta, L, h, Hb, U, fun, ksearch)
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    gainu = 10 ^ ((U.Alpha(ksearch) * log10(du2u) + U.Beta(ksearch)) / 10) / U.Noise;

    
    f = fun(gainu, gainb);
end