%% solveTwoHelixKnSBT.m
% Slender-body discretization for TWO helices offset by ±d/2 in x,
% and an optional phase shift 'delta' for the second helix.

function [xm1, ym1, zm1, xm2, ym2, zm2, ...
          T1, T2, F1, F2, ...
          vspx_xz, vspz_xz, vspx_xy, vspy_xy, vspz_xy] = ...
          solveTwoHelixKnSBT(N, a, d, R, lambda, nTurns, ...
                             Omega, V, delta, ...
                             xsp, ysp, zsp, zPlane, t)
    % 1) Discretize each helix
    phiMax  = 2*pi*nTurns;
    dphi    = phiMax / N;
    phiVec  = linspace(dphi/2, phiMax - dphi/2, N);
    deltaVal = a * exp(0.5) / 2;
    totalNodes = 2*N;
    Xall       = zeros(3, totalNodes);

    % Helix 1
    for i = 1:N
        phi = phiVec(i);
        Xall(:, i) = [
            -d/2 + R*cos(phi + Omega*t);  % x
            -R    * sin(phi + Omega*t);   % y (LH sign)
             (lambda/(2*pi)) * phi       % z
        ];
    end

    % Helix 2
    for i = 1:N
        phi = phiVec(i);
        Xall(:, N + i) = [
             d/2 + R*cos(phi + delta + Omega*t);
            -R    * sin(phi + delta + Omega*t);
             (lambda/(2*pi))*phi
        ];
    end

    % Approximate element length ds
    ds = sqrt((R*dphi)^2 + ((lambda/(2*pi))*dphi)^2);

    % 2) Tangents
    denom = sqrt(R^2 + (lambda/(2*pi))^2);
    tHat  = zeros(3, totalNodes);

    % Helix 1 tangents
    for i = 1:N
        phi = phiVec(i);
        dx  = -R*sin(phi + Omega*t);
        dy  = -R*cos(phi + Omega*t);  % LH
        dz  =  (lambda/(2*pi));
        tHat(:, i) = [dx; dy; dz] / denom;
    end

    % Helix 2 tangents
    for i = 1:N
        phi = phiVec(i);
        dx  = -R*sin(phi + delta + Omega*t);
        dy  = -R*cos(phi + delta + Omega*t); % LH
        dz  =  (lambda/(2*pi));
        tHat(:, N + i) = [dx; dy; dz] / denom;
    end

    % 3) Build system M*f = U
    M = zeros(3*totalNodes, 3*totalNodes);
    U = zeros(3*totalNodes, 1);
    I3 = eye(3);

    for n = 1:totalNodes
        rowRange = 3*(n-1) + (1:3);

        % (a) Rigid velocity
        if n <= N
            axisCenter = [-d/2; 0; 0];  % Helix 1
        else
            axisCenter = [ d/2; 0; 0];  % Helix 2
        end
        rn      = Xall(:, n);
        omegaVec= [0; 0; -Omega];  % LH rotation
        vel_n   = cross(omegaVec, rn - axisCenter) + [0; 0; V];
        U(rowRange) = vel_n;

        % (b) Local diagonal block
        tn        = tHat(:, n);
        logFactor = log((ds/2) / deltaVal);
        Kn        = logFactor * (I3 + tn*tn.');
        localDiag = (I3 - tn*tn.' + Kn)/(4*pi);
        M(rowRange, rowRange) = localDiag;

        % (c) Off-diagonal stokeslet
        for m = 1:totalNodes
            if m == n, continue; end
            colRange = 3*(m-1) + (1:3);

            rm   = Xall(:, m);
            r_nm = rn - rm;
            dist = norm(r_nm);

            S_ij = I3/dist + (r_nm*r_nm.') / (dist^3);
            M(rowRange, colRange) = M(rowRange, colRange) + ...
                                    ds * S_ij/(8*pi);
        end
    end

    % 4) Solve system
    fVec = M \ U;
    fAll = reshape(fVec, [3, totalNodes]);

    %% Compute forces, torques
    % Summation approach
    F1 = sum(fAll(3, 1:N))     * ds;
    F2 = sum(fAll(3, N+1:end)) * ds;

    % Torques
    Tvec1 = [0; 0; 0];
    Tvec2 = [0; 0; 0];
    for n = 1:2*N
        fn = fAll(:, n);
        rn = Xall(:, n);
        if n <= N
            axisCenter = [-d/2; 0; 0];
            rLocal     = rn - axisCenter;
            Tvec1      = Tvec1 + cross(rLocal, fn);
        else
            axisCenter = [ d/2; 0; 0];
            rLocal     = rn - axisCenter;
            Tvec2      = Tvec2 + cross(rLocal, fn);
        end
    end
    T1 = Tvec1(3) * ds ;
    T2 = Tvec2(3) * ds ;

    %% Return node coordinates for convenience
    xm1 = Xall(1, 1:N);  ym1 = Xall(2, 1:N);  zm1 = Xall(3, 1:N);
    xm2 = Xall(1, N+1:end); ym2 = Xall(2, N+1:end); zm2 = Xall(3, N+1:end);

    % 5) Compute velocities on x-z and x-y planes

    % Perpendicular part
    fPerpAll = zeros(3, totalNodes);
    for n = 1:totalNodes
        tn          = tHat(:, n);
        f_n         = fAll(:, n);
        fPerpAll(:, n) = (I3 - tn*tn.') * f_n;
    end

    % Evaluate velocity fields
    vspx_xz = zeros(length(xsp), length(zsp));
    vspz_xz = zeros(length(xsp), length(zsp));

    vspx_xy = zeros(length(xsp), length(ysp));
    vspy_xy = zeros(length(xsp), length(ysp));
    vspz_xy = zeros(length(xsp), length(ysp));

    for i = 1:length(xsp)
        for k = 1:length(zsp)
            [vx_val, ~, vz_val] = ...
                calculateVelocity(xsp(i), 0, zsp(k), Xall, fAll, fPerpAll, totalNodes, ds,a);
            vspx_xz(i, k) = vx_val;
            vspz_xz(i, k) = vz_val;
        end
    end

    for i = 1:length(xsp)
        for j = 1:length(ysp)
            [vx_val, vy_val, vz_val] = ...
                calculateVelocity(xsp(i), ysp(j), zPlane, Xall, fAll, fPerpAll, totalNodes, ds,a);
            vspx_xy(i, j) = vx_val;
            vspy_xy(i, j) = vy_val;
            vspz_xy(i, j) = vz_val;
        end
    end
end
