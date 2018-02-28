% *************************************************************************
%
% vortexPanel(xPts, yPts, V_inf, alpha) returns the sectional coefficient
% of lift given the vortex panel boundary positions on an airfoil 
% (specified by xPts and yPts), the free stream velocity, and the angle of 
% attack of the airfoil. 
%
%
% Syntax
% c_l = vortexPanel(xB, yB, V_inf, alpha) 
% [c_l, cpPlot] =  vortexPanel(xB, yB, V_inf, alpha) 
% 
%   Inputs:
%       - xB:      x coordinates of the panel boundaries on the airfoil.
%       - yB:      y coordinates of the panel boundaries on the airfoil. 
%       - V_inf:   Free stream velocity of uniform flow. Note that this has
%                  no effect on the results of this function. 
%       - alpha:   Angle of attack of the airfoil in the flow. 
%
%   Outputs:
%       - c_l:     sectional coefficient of lift of the airfoil in the 
%                  given flow.
%       - cpPlot:  A matrix with control point x values in the first column
%                  and their corresponding coefficients of pressure in the
%                  second
%  
%    Dependancies:
%       - In order to get the x and y coordinates of the vortex panel
%         boundaries, it is recommended that NACAairfoilPlot.m is used. 
%
%
% Created: 10/10/17 - Connor Ott
% Last Modified: 10/11/17 - Connor Ott
%
% *************************************************************************
function varargout = vortexPanel(xB, yB, V_inf, alpha)

[~, col] = size(xB);
numPanels = col - 1;
numBounds = col; 
alpha = alpha * pi/180; % [rad] - angle of attack

%% Describing Panel geometry and control point locations
xC = 0.5 * (xB(1:end - 1) + xB(2:end)); % control point locations
yC = 0.5 * (yB(1:end - 1) + yB(2:end)); 
S_panel =  sqrt((diff(xB).^2 + diff(yB).^2)); % length of each panel
theta = atan2(diff(yB), diff(xB)); % Angle of each panel WRT to horizontal

% Used frequently in equations: 
sinTheta = sin(theta);
cosTheta = cos(theta);
RHS = sin(theta - alpha)'; % right hand side of eq. 5.47 (Kuethe and Chow)
RHS(numBounds) = 0;

%% Generating systems of equations. 
% The C matrices
[C_n1, C_n2, C_t1, C_t2] = deal(zeros(numPanels, numPanels));
for i = 1:numPanels
    for j = 1:numPanels
        A = -(xC(i) - xB(j))*cosTheta(j) - (yC(i) - yB(j))*sinTheta(j);
        B = (xC(i) - xB(j))^2 + (yC(i) - yB(j))^2;
        C = sin(theta(i) - theta(j));
        D = cos(theta(i) - theta(j));
        E = (xC(i) - xB(j))*sinTheta(j) - (yC(i) - yB(j)) * cosTheta(j);
        F = log(1 + (S_panel(j)^2 + 2*A*S_panel(j))/B);
        G = atan2(E*S_panel(j), B + A*S_panel(j));
        P = (xC(i) - xB(j))*sin(theta(i) - 2*theta(j)) + ...
            (yC(i) - yB(j))*cos(theta(i) - 2*theta(j));
        Q = (xC(i) - xB(j))*cos(theta(i) - 2*theta(j)) - ...
            (yC(i) - yB(j))*sin(theta(i) - 2*theta(j));
        
        % The C matrices (to later be rewriten to make A matrices)
        C_n2(i, j) = D + 0.5*(Q*F)/S_panel(j) - (A*C + D*E)*G/S_panel(j);
        C_n1(i, j) = 0.5*D*F + C*G - C_n2(i, j);
        
        C_t2(i, j) = C + 0.5*(P*F)/S_panel(j) + (A*D - C*E)*G/S_panel(j);
        C_t1(i, j) = 0.5*C*F - D*G - C_t2(i, j);
    end
end

% Fixing up the C matrices for cases where i = j (diagonal)
repDiagC_n1 = diag(ones(1, numPanels))*(-1); % replacement diagonal matrix
repDiagC_n2 = diag(ones(1, numPanels));
C_n1 = C_n1 - diag(diag(C_n1)) + repDiagC_n1; % replacing diagonal of C_n1
C_n2 = C_n2 - diag(diag(C_n2)) + repDiagC_n2;

repDiagC_t = diag(ones(1, numPanels))*(pi/2); % same replacement for 1 & 2
C_t1 = C_t1 - diag(diag(C_t1)) + repDiagC_t; % replacing diagonal of C_t1
C_t2 = C_t2 - diag(diag(C_t2)) + repDiagC_t;

%% Condensing the C matrices into A matrices.
% A_n matrix
A_n = zeros(numBounds, numBounds);
A_n(numBounds, 1) = 1;
A_n(numBounds, numBounds) = 1;
A_n(1:numPanels, 1) = C_n1(:, 1);
A_n(1:numPanels, numBounds) = C_n2(:, numPanels);
% Combining innards of C matrices as specified in Kuethe and Chow
C_nComb = C_n1(:, 2:numPanels) + C_n2(:, 1:numPanels - 1);
A_n(1:numPanels, 2:numPanels) = C_nComb; % voila

% A_t matrix
A_t = zeros(numPanels, numBounds);
A_t(:, 1) = C_t1(:, 1);
A_t(:, numBounds) = C_t2(:, numPanels);
% Combinign C_t matrices
C_tComb = C_t1(:, 2:numPanels) + C_t2(:, 1:numPanels - 1);
A_t(1:numPanels, 2:numPanels) = C_tComb;

%% Solving the system of equations for Gamma at each control point. 
gama = A_n\(RHS); 
%% Determining Velocity using gamma

Vel_cont = zeros(1, numPanels);
for i = 1:numPanels
    Vel_cont(i) = cos(theta(i) - alpha);
    for j = 1:numBounds
        Vel_cont(i) = Vel_cont(i) + A_t(i, j)*gama(j);
    end
end

% Coefficient of pressure (Vel_cont is dimensionless)
c_pVec = 1 - (Vel_cont).^2;

%% Calculating circulation from velocity and panel length
GamTot = sum(Vel_cont.*S_panel);
c = xB(1);
varargout{1} = 2*GamTot/c;

if nargout == 2
    varargout{2} = [xC', c_pVec'];
end

end












