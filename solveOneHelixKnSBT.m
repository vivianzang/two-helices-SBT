%% solveOneHelixKnSBT_forceTorque.m
% Slender-body discretization for a SINGLE helix.

function [Fz, Tz, vspx_xz, vspz_xz, vspz_xy] = ...
    solveOneHelixKnSBT(N, R, a, lambda, nTurns, Omega, V, ...
                       xsp, ysp, zsp, zPlane, t)

    %% 1) Discretization
    phiMax  = 2*pi * nTurns;
    dphi    = phiMax / N;
    phiVec  = linspace(dphi/2, phiMax - dphi/2, N);
    deltaVal = a * exp(0.5) / 2;
    Xall    = zeros(3, N);
    for i = 1:N
        phi = phiVec(i);
        Xall(1, i) = R * cos(phi + Omega*t);
        Xall(2, i) = -R * sin(phi + Omega*t);   % LH
        Xall(3, i) = (lambda / (2*pi)) * phi;
    end

    % approximate ds
    ds = sqrt((R*dphi)^2 + ((lambda/(2*pi))*dphi)^2);

    %% 2) Local tangents
    denom = sqrt(R^2 + (lambda/(2*pi))^2);
    tHat  = zeros(3, N);
    for i = 1:N
        phi  = phiVec(i);
        dx   = -R * sin(phi + Omega*t);
        dy   = -R * cos(phi + Omega*t);  % LH
        dz   =  (lambda / (2*pi));
        tHat(:, i) = [dx; dy; dz] / denom;
    end

    %% 3) Build M f = U
    totalNodes = N;
    M = zeros(3*N, 3*N);
    U = zeros(3*N, 1);
    I3 = eye(3);

    for n = 1:totalNodes
        rowRange = 3*(n-1) + (1:3);

        % (a) boundary velocity
        rn      = Xall(:, n);
        omegaVec= [0; 0; -Omega];  % LH
        vel_n   = cross(omegaVec, rn) + [0; 0; V];
        U(rowRange) = vel_n;

        % (b) local diagonal
        tn        = tHat(:, n);
        logFactor = log(ds / (2*deltaVal));
        Kn        = logFactor * (I3 + tn*tn.');
        localDiag = (I3 - tn*tn.' + Kn)/(4*pi);
        M(rowRange, rowRange) = localDiag;

        % (c) off-diagonal stokeslet
        for m = 1:totalNodes
            if m == n, continue; end
            colRange = 3*(m-1) + (1:3);

            rm   = Xall(:, m);
            r_nm = rn - rm;
            dist = norm(r_nm);

            S_ij = I3/dist + (r_nm*r_nm.')/(dist^3);
            M(rowRange, colRange) = M(rowRange, colRange) + ...
                                    ds*S_ij/(8*pi);
            
        end
    end

    %% 4) Solve Matrix
    fVec = M \ U;
    fAll = reshape(fVec, [3, N]);

    %% 5) Compute total force & torque
    % Force
    Fnet = sum(fAll, 2); 
    Fz   = Fnet(3) * ds;

    % Torque (about z)
    Tvec = [0; 0; 0];
    for n = 1:N
        rn = Xall(:, n);
        fn = fAll(:, n);
        Tvec = Tvec + cross(rn, fn);
    end
    Tz = Tvec(3) * ds;

    %% 6) For velocity field
    fPerpAll = zeros(3, N);
    for n = 1:N
        tn = tHat(:, n);
        fn = fAll(:, n);
        fPerpAll(:, n) = (I3 - tn*tn.') * fn;  % perpendicular
    end

    % Evaluate velocity fields
    vspx_xz = zeros(length(xsp), length(zsp));
    vspz_xz = zeros(length(xsp), length(zsp));
    vspz_xy = zeros(length(xsp), length(ysp));

    % x-z plane (y=0)
    for i = 1:length(xsp)
        for k = 1:length(zsp)
            [vx_val, ~, vz_val] = ...
                calculateVelocity(xsp(i), 0, zsp(k), Xall, fAll, fPerpAll, totalNodes, ds,a);
            vspx_xz(i, k) = vx_val;
            vspz_xz(i, k) = vz_val;
        end
    end

    % x-y plane at z = zPlane
    for i = 1:length(xsp)
        for j = 1:length(ysp)
            [~, ~, vz_val] = ...
                calculateVelocity(xsp(i), ysp(j), zPlane, Xall, fAll, fPerpAll, totalNodes, ds,a);
            vspz_xy(i, j) = vz_val;
        end
    end
end
