function Rhos0 = findCriticalPoints2DK(Rhos, ks, PosUAV, PosUE, PosBS, U, fun)
% Version 2DK, search in 2D case, not to optimize UAV height
% Version K, for K segment case

    L = norm(PosUE(1:2) - PosBS(1:2), 2);
    rho = norm(PosUAV(1:2) - PosUE(1:2), 2);
    Uvec = (PosBS(1:2) - PosUE(1:2)) / norm(PosBS(1:2) - PosUE(1:2), 2);
    Posdelta = PosUAV(1:2) - PosUE(1:2);

    gamma= Posdelta(:).' * Uvec(:) / rho;
    theta = sign(Posdelta(2) * Uvec(1) - Posdelta(1) * Uvec(2)) * real(acos(gamma));

    h = PosUAV(3);
    Hb = PosBS(3);

    if abs(theta) > 1e-6
        error('theta should be zero!');
    end

    Rhos0 = zeros(1, U.K);
    
    for k = 1:U.K
        I = find(ks == k);
        rhomin = min(Rhos(I));
        rhomax = max(Rhos(I));
        if isempty(I)
            Rhos0(k) = Rhos0(k - 1);
            continue
        end
        if k == 1
            rhomin = 0;
        elseif k == U.K
            rhomax = L;
        end 
        los = k; 
            
        delta = 1e-4;
        tol = L * 1e-7;
        cnt = 0;
        MAXLOOP = 1000;
        while rhomax - rhomin > tol && cnt < MAXLOOP
            cnt = cnt + 1 ;
            rho = (rhomax + rhomin) / 2;

            F_plus = F(rho + delta, theta, L, h, Hb, U, fun, los);
            F_minus = F(rho - delta, theta, L, h, Hb, U, fun, los);
            dF_drho = (F_plus - F_minus) / (2 * delta);

            if dF_drho >= 0
                rhomax = rho;
            else
                rhomin = rho;
            end
        end

        if cnt >= MAXLOOP
            error('rho did not converge!');
            rho0 = 0;
        else
            rho0 = (rhomax + rhomin) / 2;
        end   
        Rhos0(k) = rho0;
        
    end

end

function f = F(rho, theta, L, h, Hb, U, fun, los)
    % compute the function F(\rho,\theta) in the paper

    lopt = sqrt(rho^2 + h^2);
    f = fl(lopt, rho, theta, L, h, Hb, U, fun, los);
    
end

function f = fl(l, rho, theta, L, h, Hb, U, fun, los)
    k = los;
    Db = L^2 * sin(theta)^2 + (L * cos(theta) - rho * l / sqrt(rho^2 + h^2))^2 ...
        + (l * h / sqrt(rho^2 + h^2) - Hb)^2;
    du2b = sqrt(Db);
    gainb = 10 ^ ((U.A0 * log10(du2b) + U.B0) / 10) / U.Noise;
    
    du2u = l;
    gainu = 10 ^ ((U.Alpha(k) * log10(du2u) + U.Beta(k)) / 10) / U.Noise;
    
    f = fun(gainu, gainb);
end
