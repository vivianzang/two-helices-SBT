%% timeAvgAxialVelDouble.m
% Compute time-averaged axial velocity for a double-helix configuration 
% over one rotation (0 to 2*pi).

function [vzTimeAvg, F1avg, F2avg, T1avg, T2avg] = timeAvgAxialVelDouble( ...
    N, a, d, R, lambda, nTurns, Omega, V, delta, ...
    xsp, ysp, zsp, zPlane)

    % Number of time steps for averaging
    nTimeSteps = 25;
    timeVec    = linspace(0, 2*pi, nTimeSteps);

    % Preallocate
    vzIns   = zeros(length(xsp), length(ysp), nTimeSteps);
    thrust1 = zeros(1, nTimeSteps);
    thrust2 = zeros(1, nTimeSteps);
    torque1 = zeros(1, nTimeSteps);
    torque2 = zeros(1, nTimeSteps);

    for i = 1:nTimeSteps
        t = timeVec(i);
        
        [~, ~, ~, ~, ~, ~, T1, T2, F1, F2, ~, ~, ~, ~, vspz2d] = ...
            solveTwoHelixKnSBT(N, a, d, R, lambda, nTurns, ...
                               Omega, V, delta, ...
                               xsp, ysp, zsp, zPlane, t);
        
        vzIns(:, :, i) = vspz2d;
        thrust1(i)     = F1;
        thrust2(i)     = F2;
        torque1(i)     = T1;
        torque2(i)     = T2;
    end

    % Time average
    vzTimeAvg = mean(vzIns, 3);
    F1avg     = mean(thrust1);
    F2avg     = mean(thrust2);
    T1avg     = mean(torque1);
    T2avg     = mean(torque2);
end
