%% timeAvgAxialVzSingle.m
% Compute time-averaged axial velocity for a single helix configuration.

function [Fz, Tz, vzTimeAvg] = timeAvgAxialVzSingle( ...
    N, R, a, lambda, nTurns, Omega, V, ...
    xsp, ysp, zsp, zPlane)

    % Number of time steps for averaging
    nTimeSteps = 25;
    timeVec    = linspace(0, 2*pi, nTimeSteps);

    % Preallocate
    vzIns = zeros(length(xsp), length(ysp), nTimeSteps);
    FzArr = zeros(1, nTimeSteps);
    TzArr = zeros(1, nTimeSteps);

    for i = 1:nTimeSteps
        t = timeVec(i);
        
        [Fz_i, Tz_i, vspx_xz, vspz_xz, vspz_xy] = ...
            solveOneHelixKnSBT(...
                N, R, a, lambda, nTurns, Omega, V, ...
                xsp, ysp, zsp, zPlane, t);

        vzIns(:, :, i) = vspz_xy;
        FzArr(i)       = Fz_i;
        TzArr(i)       = Tz_i;
    end

    % Time averages
    vzTimeAvg = mean(vzIns, 3);
    Fz        = mean(FzArr);
    Tz        = mean(TzArr);
end
