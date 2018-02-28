% *************************************************************************
%
% plotAirfoilFlow(c, alpha, V_inf, P_inf, rho_inf, N)
%
% Plots the pressure contours and stream lines around a symmetric airfoil
% using superpostion of N vortices given the dimensions of the airfoil,
% flow properties, and number of vortices used for calculation.
%
%   Inputs:
%       c - Chord Length [m]
%       alpha - Angle of attack [deg]
%       V_inf - Free stream velocity [m/s]
%       P_inf - Free stream pressure [Pa]
%       rho_inf - Free stream air density [kg/m^3]
%       N - Number of vortices used for model
%   Outputs:
%       Plot of airfoil, streamlines, and pressure contour lines.
%   Dependencies:
%       None
%
% Created: 9/28/17 - Connor Ott
% Last Modified: 10/06/17 - Connor Ott
%
% *************************************************************************

function[] = plotAirfoilFlow(c, alpha, V_inf, P_inf, rho_inf, N)

tic
alpha = alpha * pi/180; % [Rad] - Angle of attack

% Starting points (limited to just the airfoil now)
xmin = -c/2;
xmax = 1.5*c;
ymin = -c;
ymax = c;
numMesh = 100;
numMesh_vis = 16;

% Defining mesh for flow field
[xMat,yMat]=meshgrid(linspace(xmin,xmax,numMesh),...
            linspace(ymin,ymax,numMesh));
[xMat_vis, yMat_vis] = meshgrid(linspace(xmin,xmax,numMesh_vis),...
                       linspace(ymin,ymax,numMesh_vis));
%% Defining Stream Functions and Velocity Potentials

radius = @(x, y, x_pos) ((x-x_pos).^2 + y.^2).^(1/2); % Radius in cartesean

Z = alpha * V_inf * c / (N * pi); % Consolidating Constants
x_posVec = linspace(0, c, N+1);

% Uniform Flow
psi_U = V_inf * (cos(alpha)* yMat - sin(alpha) * xMat);
phi_U = V_inf * (cos(alpha)* xMat + sin(alpha) * yMat);
u_U = V_inf * cos(alpha); % {u} component of {U}niform flow velocity
v_U = V_inf * sin(alpha); 

% Initial uniform flow field properties (to be added to vortices)
psiTotal = psi_U;
phiTotal = phi_U;
u_tot = u_U;
v_tot = v_U;

u_totVis = u_U;
v_totVis = v_U;

wait = waitbar(N, sprintf('Calculating Fields Using %.e Vortices', N));
set(findall(wait,'type','text'),'Interpreter','none', ...
                                'FontName', 'Segoe UI' ,...
                                'FontSize', 9');
wait.Name = 'MATLAB';

% Summing vortices for streamlines and potential
for i = 1:N
    x_i = x_posVec(i+1); 
    
    % Vortex Stream Functions
    psi_VortTemp = Z * sqrt((1 - x_i/c)/(x_i/c)) .* ...
                   log(radius(xMat, yMat, x_i));
    psiTotal = psiTotal + psi_VortTemp;
    
    % Vortex Velocity Potential
    phi_VortTemp = - Z * sqrt((1 - x_i/c)/(x_i/c)) .* atan2(yMat, (xMat-x_i));
    phiTotal = phiTotal + phi_VortTemp;
    
    % Vortex Velocity Field
    u_Vtemp =  Z * yMat * sqrt((1 - x_i/c)/(x_i/c)) .* ...
        radius(xMat, yMat, x_i).^(-2);
    v_Vtemp = - Z * (xMat - x_i) * sqrt((1 - x_i/c)/(x_i/c)) .*  ...
        radius(xMat, yMat, x_i).^(-2);
    u_tot = u_tot + u_Vtemp;
    v_tot = v_tot + v_Vtemp;
    
    % Flow Velocity visualization field
    u_VtVis = Z * yMat_vis * sqrt((1 - x_i/c)/(x_i/c)) .* ...
        radius(xMat_vis, yMat_vis, x_i).^(-2);
    v_VtVis = - Z * (xMat_vis - x_i) * sqrt((1 - x_i/c)/(x_i/c)) .*  ...
        radius(xMat_vis, yMat_vis, x_i).^(-2);
    u_totVis = u_totVis + u_VtVis;
    v_totVis = v_totVis + v_VtVis;
    
    waitbar(i/N)
end
close(wait)

%% Velocity Field -> Pressure field
VelTot_comp = {u_tot, v_tot};
VelTot_mag = sqrt(VelTot_comp{1}.^2 + VelTot_comp{2}.^2);

% Bernoulli's
P_field = P_inf + 0.5*rho_inf*(V_inf^2 - VelTot_mag.^2);

% Sectional Pressure Coefficient
Q_inf = 0.5 * rho_inf * V_inf^2;
c_p = (P_field - P_inf)./(0.5 * rho_inf * V_inf^2 * c);

%% Determine color levels for contours - Adapted from Lifting_Cylinder.m
% defines the color levels -> trial and error to find a good representation
levminPsi = psiTotal(1,numMesh); 
levmaxPsi = psiTotal(numMesh,numMesh/2);
levelsPsi = linspace(levminPsi,levmaxPsi, 40)';

levminPhi = min(min(phiTotal)); % defines the color levels
levmaxPhi = max(max(phiTotal));
levelsPhi = linspace(levminPhi,levmaxPhi,40)';

levminPres = min(min(P_inf-2.75*Q_inf)); % (Thanks Marshall)
levmaxPres = max(max(P_inf+1.5*Q_inf));
levelsPres = linspace(levminPres, levmaxPres, 20)';



%% Plotting it all
set(0, 'defaulttextInterpreter', 'latex')

% Equipotential and Streamlines *******************************************
figure
hold on

contour(xMat, yMat, psiTotal, levelsPsi)
plot([0 c], [0 0], 'k', 'linewidth', 3) % the "airfoil"
title(sprintf(['Stream Lines at $%.f',...
               '^{\\circ}$ Angle of Attack'], alpha * 180/pi));
xlabel('x Position [m]');
ylabel('y Position [m]');

set(gca, 'TickLabelInterpreter', 'latex')
coBar1 = colorbar;
coBar1.Label.String = 'Stream Line Values [~]';
coBar1.Label.Interpreter = 'latex';
coBar1.Label.FontSize = 12;
coBar1.TickLabelInterpreter = 'latex';
axis equal
axis([xmin xmax ymin/2 ymax/2])
hold off
saveas(gcf,'stream.png')

% Equipotential Lines *****************************************************
figure
hold on

contour(xMat, yMat, phiTotal, levelsPhi)
plot([0 c], [0 0], 'k', 'linewidth', 3) % the "airfoil"
title(sprintf(['Equipotential Lines at $%.f',...
               '^{\\circ}$ Angle of Attack'], alpha * 180/pi));
xlabel('x Position [m]');
ylabel('y Position [m]');

set(gca, 'TickLabelInterpreter', 'latex')
coBar1 = colorbar;
coBar1.Label.String = 'Velocity Potential [~]';
coBar1.Label.Interpreter = 'latex';
coBar1.Label.FontSize = 12;
coBar1.TickLabelInterpreter = 'latex';
axis equal
axis([xmin xmax ymin/2 ymax/2])
hold off
saveas(gcf, 'equipot.png')

% Pressure Contour Lines **************************************************
figure
hold on

contourf(xMat, yMat, P_field, levelsPres)
colormap jet
plot([0 c], [0 0], 'k', 'linewidth', 3) % the "airfoil"

title(sprintf(['Pressure Field Contours at $%.f ^{\\circ}$ Angle of',...
               'Attack'], alpha * 180/pi));
xlabel('x Position [m]');
ylabel('y Position [m]');

coBar2 = colorbar;
coBar2.Label.String = 'Pressure [Pa]';
coBar2.Label.FontSize = 12;
coBar2.Label.Interpreter = 'latex';
coBar2.TickLabelInterpreter = 'latex';
set(gca, 'TickLabelInterpreter', 'latex')
axis equal
axis([xmin xmax ymin/2 ymax/2])
hold off
saveas(gcf,'Pcontour.png')


% Velocity Field graph (Low Res) ******************************************
figure
hold on

[~, l] = contourf(xMat, yMat, phiTotal, levelsPhi);
set(l,'LineColor','none') % no isolines
quiver(xMat_vis,yMat_vis,u_totVis,v_totVis, 'k')
plot([0 c], [0 0], 'k', 'linewidth', 3) % the "airfoil"

title(sprintf('Velocity Field at $%.f ^{\\circ}$ Angle of Attack', ...
    alpha * 180/pi));
xlabel('x Position [m]');
ylabel('y Position [m]');
axis equal
axis([xmin xmax ymin/2 ymax/2])

coBar = colorbar;
coBar.Label.String = 'Velocity Potential [~]';
coBar.Label.FontSize = 12;
coBar.Label.Interpreter = 'latex';
coBar.TickLabelInterpreter = 'latex';
hold off

saveas(gcf,'V_field.png')
% *************************************************************************
runTime = toc;
if runTime > 20
    fprintf(['Well if you didn''t want to wait so long why''d you use',...
            'such a large N?\n']);
    fprintf('Elapsed Time is %.1f', runTime);
end



end





