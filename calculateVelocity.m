%% calculateVelocity.m
% Computes the velocity at a point (x, y, z) due to a set of forces 'f1'
% and a set of perpendicular (dipole) terms 'f2' along a curve X.

function [vx, vy, vz] = calculateVelocity( ...
    x, y, z, X, f1, f2, nRes, ds,a)

    vi      = [0; 0; 0];
    xflow   = [x; y; z];

    for dp = 1:nRes
        rDiff  = xflow - X(:, dp);
        dist   = norm(rDiff);

        % Stokeslet
        S_ij = eye(3)/dist + (rDiff*rDiff')/(dist^3);

        % Dipole
        D_ij = -eye(3)/dist^3 + 3*(rDiff*rDiff')/(dist^5);

        % Combine
        integrand = (S_ij * f1(:, dp))/(8*pi) - (D_ij * a^2 * f2(:, dp)) / (16*pi);

        % Accumulate
        vi = vi + integrand * ds;
    end

    vx = vi(1);
    vy = vi(2);
    vz = vi(3);
end
