% *************************************************************************
%
% NOplotAirfoilFlow(c, alpha, V_inf, P_inf, rho_inf, N)
%
% Gives pressure and velocity fields only for error analysis. For plots,
% see plotAirfoilFlow.
%
%   Inputs:
%       c - Chord Length [m]
%       alpha - Angle of attack [deg]
%       V_inf - Free stream velocity [m/s]
%       P_inf - Free stream pressure [Pa]
%       rho_inf - Free stream air density [kg/m^3]
%       N - Number of vortices used for model
%   Outputs: (Dependant on nargout)
%       P_field - Pressure distribution
%       VelTot_mag - Velocity magnitude field
%       xMat - Mesh of x values
%       yMat - Mesh of y values
%   Dependencies:
%       None
%
% Created: 9/28/17 - Connor Ott
% Last Modified: 10/06/17 - Connor Ott
%
% *************************************************************************

function[varargout] = NOplotAirfoilFlow(c, alpha, V_inf, P_inf, rho_inf, N)


alpha = alpha * pi/180; % [Rad] - Angle of attack

% Starting points (limited to just the airfoil now)
xmin = -c/2;
xmax = 1.5*c;
ymin = -c;
ymax = c;
numMesh = 100;

% Defining mesh for flow field
[xMat,yMat] = meshgrid(linspace(xmin,xmax,numMesh),...
              linspace(ymin,ymax,numMesh));
%% Defining Stream Functions and Velocity Potentials

radius = @(x, y, x_pos) ((x-x_pos).^2 + y.^2).^(1/2); % Radius in cartesean

Z = alpha * V_inf * c / (N * pi); % Consolidating Constants
x_posVec = linspace(0, c, N+1);

% Uniform Flow
psi_U = V_inf * (cos(alpha)* yMat - sin(alpha) * xMat);
phi_U = V_inf * (cos(alpha)* xMat + sin(alpha) * yMat);
u_U = V_inf * cos(alpha); % u component of Uniform flow velocity
v_U = V_inf * sin(alpha);

% Initial uniform flow field properties (to be added to vortices)
psiTotal = psi_U;
phiTotal = phi_U;
u_tot = u_U;
v_tot = v_U;

% Summing vortices for streamlines and potential
wait = waitbar(N, sprintf('Calculating Fields Using %.1e Vortices', N));
set(findall(wait,'type','text'),'Interpreter','none', ...
                                'FontName', 'Segoe UI' ,...
                                'FontSize', 9');
wait.Name = 'MATLAB';

for i = 1:N
    x_i = x_posVec(i+1); %
    
    % Vortex Stream Functions
    psi_VortTemp = Z * sqrt((1 - x_i/c)/(x_i/c)) .* ...
                   log(radius(xMat, yMat, x_i));
    psiTotal = psiTotal + psi_VortTemp;
    
    % Vortex Velocity Potential
    phi_VortTemp = - Z * sqrt((1 - x_i/c)/(x_i/c)) .* ...
                   atan2(yMat, (xMat-x_i));
    phiTotal = phiTotal + phi_VortTemp;
    
    % Vortex Velocity Field
    u_Vtemp =  Z * yMat * sqrt((1 - x_i/c)/(x_i/c)) .* ...
        radius(xMat, yMat, x_i).^(-2);
    v_Vtemp = - Z * (xMat - x_i) * sqrt((1 - x_i/c)/(x_i/c)) .*  ...
        radius(xMat, yMat, x_i).^(-2);
    u_tot = u_tot + u_Vtemp;
    v_tot = v_tot + v_Vtemp;
    
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



%% Outputs in case I ask for them in errorAnalysis.m
switch nargout
    case 2
        varargout{1} = P_field;
        varargout{2} = VelTot_mag;
    case 4
        varargout{1} = P_field;
        varargout{2} = VelTot_mag;
        varargout{3} = xMat;
        varargout{4} = yMat;
    otherwise
        % Do not assign
end


end





