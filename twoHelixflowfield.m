function twoHelixflowfield(N, a, d, R, lambda, nTurns, Omega, V, delta,velocityType)
% Spatial grid definitions
xsp = linspace(-10,10,51);
ysp = linspace(-10,10,50);
zsp = linspace(-3,30,61);
[Ysp3, Xsp3] = meshgrid(ysp, xsp);  % For x-y plane plots
[Zsp2, Xsp2] = meshgrid(zsp, xsp);  % For x-z plane plots
zPlane = 15;
% Animation parameters
nFrames = 25;           % Number of frames (e.g., one rotation)
timeVals = linspace(0,2*pi,nFrames);
% Preallocate velocity arrays (for clarity)
vspx_2d = zeros(length(xsp), length(zsp), nFrames);
vspz_2d = zeros(length(xsp), length(zsp), nFrames);
vx_xy   = zeros(length(xsp), length(ysp), nFrames);
vy_xy   = zeros(length(xsp), length(ysp), nFrames);
vz_xy   = zeros(length(xsp), length(ysp), nFrames);

    % Create figure for animation
    % fig = figure;
    for i = 1:nFrames
        t = timeVals(i);
        
        % Compute the double helix flow field using your external function.
        % The function is assumed to output:
        %   xm1, ym1, zm1 : Centerline of helix 1
        %   xm2, ym2, zm2 : Centerline of helix 2
        %   vspx_2d, vspz_2d: x- and z-components (for x-z plane plots)
        %   vx_xy, vy_xy, vz_xy: x-, y-, z-components (for x-y plane plots)
        [xm1, ym1, zm1, xm2, ym2, zm2, ...
          ~, ~, ~, ~, ...
          vspx_2d(:,:,i), vspz_2d(:,:,i), vx_xy(:,:,i), vy_xy(:,:,i), vz_xy(:,:,i)] = ...
          solveTwoHelixKnSBT(N, a, d, R, lambda, nTurns, ...
                             Omega, V, delta, ...
                             xsp, ysp, zsp, zPlane, t);
        % Depending on the chosen velocity type, plot the appropriate field.
        switch lower(velocityType)
            case 'axial'
                % Plot x-z plane axial velocity using vspz_2d
                contourf(Xsp2, Zsp2, vspz_2d(:,:,i), 200, ...
                    'FaceAlpha', 0.9, 'LineColor', 'none');
                xlabel('x/R','Interpreter','latex');
                ylabel('z/R','Interpreter','latex');
                axis equal;
                xlim([-10 10]); ylim([0 20]);
                hold on;
                plot(xm1, zm1, 'k--');
                plot(xm2, zm2, 'k--');
                hold off;
                colorbar; 
                caxis([-0.4, 0.2]);
                
            case 'transverse'
                % Plot x-z plane transverse velocity using vspx_2d
                contourf(Xsp2, Zsp2, vspx_2d(:,:,i), 200, ...
                    'FaceAlpha', 0.9, 'LineColor', 'none');
                xlabel('x/R','Interpreter','latex');
                ylabel('z/R','Interpreter','latex');
                axis equal;
                xlim([-10 10]); ylim([0 20]);
                hold on;
                plot(xm1, zm1, 'k--');
                plot(xm2, zm2, 'k--');
                hold off;
                colorbar; 
                caxis([-0.5, 0.5]);
                
            case 'azimuthal'
                % Plot x-y plane azimuthal velocity magnitude
                contourf(Xsp3, Ysp3, sqrt(vx_xy(:,:,i).^2 + vy_xy(:,:,i).^2), 100, 'FaceAlpha', 0.9, 'LineStyle', "none");
                axis equal;
                xlim([-8 8]); ylim([-5 5]);
                xlabel('x/R','Interpreter','latex');
                ylabel('y/R','Interpreter','latex');
                colorbar; 
                caxis([0 0.8]);
                l = streamslice(Xsp3', Ysp3', vx_xy(:,:,i)', vy_xy(:,:,i)');
                set(l, 'Color', 'w');
                hold on;
                plot(xm1(1:20), ym1(1:20), 'k--');
                plot(xm2(1:20), ym2(1:20), 'k--');
                hold off;
                
            otherwise
                error('Invalid velocity type. Use "axial", "transverse", or "azimuthal".');
        end
        
        pause(0.1); % Pause to control the animation speed
        drawnow;
    end
end