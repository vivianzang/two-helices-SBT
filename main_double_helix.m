%% main_double_helix.m
% This is the main script that:
% 1) sets parameters 
% 2) animate flow fiel ('axial','transverse' or 'azimuthal')
% 3-5) runs the double-helix and
% 6) single-helix computations, and processes or visualizes results.

clear; close all; clc;

%% 1) Set up parameters
R       = 1; %helix radius
a       = 0.1 * R; %filament radius
L       = 36 * R; %axial length
lambda  = 12 * R; %pitch length
nTurns  = L / lambda;   
N       = 60;          % discretization points per helix
Omega   = 1;
V       = 0;
mu      = 1;
deltaVal= a * exp(0.5) / 2;
d = 5; %separation distance
delta = 0; %phase difference

%% 2) Animate double helix flow field visualization
twoHelixflowfield(N, a, d, R, lambda, nTurns, Omega, V, delta,'azimuthal')%'axial','transverse' or 'azimuthal'

%% 3) Arrays for your parametric sweeps
disArray = [2.05 4 6 8 10];  % separation distances
delArray = linspace(0, pi, 13);

%% 4) Spatial grids for velocity sampling and Initialize storage for results
aPlot     = 10;                        % half-length for x, y space
xsp    = linspace(-aPlot, aPlot, 51);
ysp    = linspace(-aPlot, aPlot, 50);
zsp    = linspace(0, 36, 37);
zPlane = 15;
% For circular masking in x-y plane
[Xsp3, Ysp3] = meshgrid(xsp, ysp);
[thetaVals, rVals] = cart2pol(Xsp3, Ysp3);

% Create a mask for radial distance > aPlot
maskR = (rVals > aPlot);

fluxValPolar = zeros(length(disArray), length(delArray));
avgThrust    = zeros(length(disArray), length(delArray));
avgTorque    = zeros(length(disArray), length(delArray));

%% 5) Double Helix: Loop over 'disArray' (d) and 'delArray' (delta)
for q = 1:length(disArray)
    for p = 1:length(delArray)
        d     = disArray(q);
        delta = delArray(p);

        [vzTimeAvg, F1avg, ~, T1avg, ~] = timeAvgAxialVelDouble( ...
            N, a, d, R, lambda, nTurns, Omega, V, delta, ...
            xsp, ysp, zsp, zPlane);

        % Mask out region where r > aPlot (circular region in the plane)
        tempVz = vzTimeAvg;
        tempVz(maskR) = 0;

        % Integrate velocity to get flux
        fluxValPolar(q, p) = trapz(ysp, trapz(xsp, tempVz, 1));

        % Record thrust / torque
        avgThrust(q, p) = F1avg;
        avgTorque(q, p) = T1avg;
    end
end

%% 6) Single Helix example
[FzSingle, TzSingle, vzTimeAvgSingle] = timeAvgAxialVzSingle( ...
    N, R, a, lambda, nTurns, Omega, V, xsp, ysp, zsp, zPlane);

% Integrate for flux with circular mask
tempVzSingle = vzTimeAvgSingle;
tempVzSingle(maskR) = 0;
fluxValSinglePolar = trapz(ysp, trapz(xsp, tempVzSingle, 1));

%% 7) Display or save results
disp('--- Double Helix Results ---');
disp('Thrust:');   disp(avgThrust* 100 * 0.5 * (0.01)^2 );
disp('Torque:');   disp(avgTorque* 100 * 0.5 * (0.01)^3 );

disp('--- Single Helix Results ---');
disp('Single Thrust:'); disp(FzSingle* 100 * 0.5 * (0.01)^2);
disp('Single Torque:'); disp(TzSingle* 100 * 0.5 * (0.01)^3);

% for i = 1:length(disArray)
%     plot(delArray,avgTorque(i,:)./TzSingle,'--')
%     plot(delArray,avgThrust(i,:)./FzSingle,'-.')
%     hold on
% end
% yline(1,'-')
% xline(pi/2,'-')
% hold off

% for i = 1:length(disArray)
%     plot(delArray,fluxValPolar(i,:)./2/fluxValSinglePolar,'--')
%     hold on
% end
% yline(1,'-')
% xline(pi/2,'-')
% hold off

% Save data to .mat files
% save('doubleHelixResults.mat', 'fluxValPolar','avgThrust','avgTorque');
% save('singleHelixResults.mat', 'FzSingle','TzSingle','fluxValSinglePolar');
