% *************************************************************************
%
% vortexPanel_flap(xB1, yB1, xB2, yB2, V_inf, alpha) 
% 
% Returns the sectional coefficient of lift given the vortex panel boundary 
% positions on an airfoil and its flap, the free stream velocity, the 
% angle of attack of the airfoil, and the deflection angle of the flap. 
%
% Syntax
% c_l = vortexPanel(xB1, yB1, xB2, yB2 V_inf, alpha) 
% [c_l, cpPlot] =  vortexPanel(xB1, yB1, xB2, yB2 V_inf, alpha) 
% 
%   Inputs:
%       xB1    - x coordinates of the panel boundaries on primary body.
%       yB1    - y coordinates of the panel boundaries on primary body. 
%       xB2    - x coordinates of panel boundaries on secondary body.
%       xB2    - y coordinates of panel boundaries on secondary body.
%       V_inf  - Free stream velocity of uniform flow.
%       alpha  - Angle of attack of the airfoil in the flow. 
%
%   Outputs:
%       c_l    - sectional coefficient of lift of the airfoil in the 
%                given flow.
%       cpPlot - A matrix with control point x values in the first column
%                and their corresponding coefficients of pressure in the
%                second
%  
%    Dependancies:
%       In order to get the x and y coordinates of the vortex panel
%       boundaries, it is recommended that NACAairfoilPlot.m is used. 
%
%
% Created: 10/10/17 - Connor Ott
% Last Modified: 10/18/17 - Connor Ott
%
% *************************************************************************
function varargout = vortexPanel_flap(xB_P, yB_P, xB_S, yB_S, V_inf, alpha)

[~, numBounds] = size(xB_P);
numPanels = numBounds - 1;
alpha = alpha * pi/180; % [rad] - angle of attack

%% Describing Panel geometry and control point location 
% Control Points for (primary and secondary)
xC_P = 0.5 * (xB_P(1:end - 1) + xB_P(2:end)); 
yC_P = 0.5 * (yB_P(1:end - 1) + yB_P(2:end)); 
xC_S = 0.5 * (xB_S(1:end - 1) + xB_S(2:end)); 
yC_S = 0.5 * (yB_S(1:end - 1) + yB_S(2:end)); 

% Length of each panel (secondary and primary)
S_panelP = sqrt((diff(xB_P).^2 + diff(yB_P).^2));
S_panelS = sqrt((diff(xB_S).^2 + diff(yB_S).^2)); 

% Angle of panel WRT to pos-x axis (secondary and primary)
thetaP = atan2(diff(yB_P), diff(xB_P)); 
thetaS = atan2(diff(yB_S), diff(xB_S)); 

% Expressions used frequently in equations
sinThetaS = sin(thetaS);
sinThetaP = sin(thetaP);
cosThetaS = cos(thetaS);
cosThetaP = cos(thetaP);

% Right hand side of eq. 5.47 (Kuethe and Chow) (primary and secondary)
RHS_P = sin(thetaP - alpha)'; 
RHS_P(numBounds) = 0;
RHS_S = sin(thetaS - alpha)';
RHS_S(numBounds) = 0;

%% Generating systems of equations (1/4) - Secondary on Primary
% The C matrices
[C_n1.SonP, C_n2.SonP, C_t1.SonP, C_t2.SonP] = ...
                                         deal(zeros(numPanels, numPanels));
                                     
% Control points on the primary, vortices on the secondary
% I - Primary
% J - Secondary
xC = xC_P; 
xB = xB_S;
yC = yC_P;
yB = yB_S;
cosTheta = cosThetaS; % J - Secondary
sinTheta = sinThetaS;
S_panel = S_panelS; 

for i = 1:numPanels
    for j = 1:numPanels
        A = -(xC(i) - xB(j))*cosTheta(j) - (yC(i) - yB(j))*sinTheta(j);
        B = (xC(i) - xB(j))^2 + (yC(i) - yB(j))^2;
        C = sin(thetaP(i) - thetaS(j));
        D = cos(thetaP(i) - thetaS(j));
        E = (xC(i) - xB(j))*sinTheta(j) - (yC(i) - yB(j)) * cosTheta(j);
        F = log(1 + (S_panel(j)^2 + 2*A*S_panel(j))/B);
        G = atan2(E*S_panel(j), B + A*S_panel(j));
        P = (xC(i) - xB(j))*sin(thetaP(i) - 2*thetaS(j)) + ...
            (yC(i) - yB(j))*cos(thetaP(i) - 2*thetaS(j));
        Q = (xC(i) - xB(j))*cos(thetaP(i) - 2*thetaS(j)) - ...
            (yC(i) - yB(j))*sin(thetaP(i) - 2*thetaS(j));
        
        % The C matrices (to later be rewriten to make A matrices)
        C_n2.SonP(i, j) = D + 0.5*(Q*F)/S_panel(j) - ...
                         (A*C + D*E)*G/S_panel(j);
        C_n1.SonP(i, j) = 0.5*D*F + C*G - C_n2.SonP(i, j);
        
        C_t2.SonP(i, j) = C + 0.5*(P*F)/S_panel(j) + ...
                         (A*D - C*E)*G/S_panel(j);
        C_t1.SonP(i, j) = 0.5*C*F - D*G - C_t2.SonP(i, j);
    end
end

% NOT FIXING C BECAUSE NO PANELS OVERLAP WITH CONTROL POINTS

% ******* Condensing the C matrices into A matrices. *******
% A_n matrix
A_n.SonP = zeros(numBounds, numBounds);
A_n.SonP(1:numPanels, 1) = C_n1.SonP(:, 1);
A_n.SonP(1:numPanels, numBounds) = C_n2.SonP(:, numPanels);
% Combining innards of C matrices as specified in Kuethe and Chow
C_nComb = C_n1.SonP(:, 2:numPanels) + C_n2.SonP(:, 1:numPanels - 1);
A_n.SonP(1:numPanels, 2:numPanels) = C_nComb; % voila

% A_t matrix
A_t.SonP = zeros(numPanels, numBounds);
A_t.SonP(:, 1) = C_t1.SonP(:, 1);
A_t.SonP(:, numBounds) = C_t2.SonP(:, numPanels);
% Combinign C_t matrices
C_tComb = C_t1.SonP(:, 2:numPanels) + C_t2.SonP(:, 1:numPanels - 1);
A_t.SonP(1:numPanels, 2:numPanels) = C_tComb;

%% Generating systems of equations (2/4) - Primary on Secondary
% The C matrices
[C_n1.PonS, C_n2.PonS, C_t1.PonS, C_t2.PonS] = ...
                                         deal(zeros(numPanels, numPanels));
                                     
% Control points on the secondary, vortices on the primary
% J - primary
% I - secondary
xC = xC_S;
xB = xB_P;
yC = yC_S;
yB = yB_P;
cosTheta = cosThetaP; 
sinTheta = sinThetaP;
S_panel = S_panelP;

for i = 1:numPanels
    for j = 1:numPanels
        A = -(xC(i) - xB(j))*cosTheta(j) - (yC(i) - yB(j))*sinTheta(j);
        B = (xC(i) - xB(j))^2 + (yC(i) - yB(j))^2;
        C = sin(thetaS(i) - thetaP(j));
        D = cos(thetaS(i) - thetaP(j));
        E = (xC(i) - xB(j))*sinTheta(j) - (yC(i) - yB(j)) * cosTheta(j);
        F = log(1 + (S_panel(j)^2 + 2*A*S_panel(j))/B);
        G = atan2(E*S_panel(j), B + A*S_panel(j));
        P = (xC(i) - xB(j))*sin(thetaS(i) - 2*thetaP(j)) + ...
            (yC(i) - yB(j))*cos(thetaS(i) - 2*thetaP(j));
        Q = (xC(i) - xB(j))*cos(thetaS(i) - 2*thetaP(j)) - ...
            (yC(i) - yB(j))*sin(thetaS(i) - 2*thetaP(j));
        
        % The C matrices (to later be rewriten to make A matrices)
        C_n2.PonS(i, j) = D + 0.5*(Q*F)/S_panel(j) - ...
                         (A*C + D*E)*G/S_panel(j);
        C_n1.PonS(i, j) = 0.5*D*F + C*G - C_n2.PonS(i, j);
        
        C_t2.PonS(i, j) = C + 0.5*(P*F)/S_panel(j) + ...
                         (A*D - C*E)*G/S_panel(j);
        C_t1.PonS(i, j) = 0.5*C*F - D*G - C_t2.PonS(i, j);
    end
end

% NOT FIXING C BECAUSE NO PANELS OVERLAP WITH CONTROL POINTS

% ******* Condensing the C matrices into A matrices. *******
% A_n matrix
A_n.PonS = zeros(numBounds, numBounds);
A_n.PonS(1:numPanels, 1) = C_n1.PonS(:, 1);
A_n.PonS(1:numPanels, numBounds) = C_n2.PonS(:, numPanels);
% Combining innards of C matrices as specified in Kuethe and Chow
C_nComb = C_n1.PonS(:, 2:numPanels) + C_n2.PonS(:, 1:numPanels - 1);
A_n.PonS(1:numPanels, 2:numPanels) = C_nComb; % voila

% A_t matrix
A_t.PonS = zeros(numPanels, numBounds);
A_t.PonS(:, 1) = C_t1.PonS(:, 1);
A_t.PonS(:, numBounds) = C_t2.PonS(:, numPanels);
% Combining C_t matrices
C_tComb = C_t1.PonS(:, 2:numPanels) + C_t2.PonS(:, 1:numPanels - 1);
A_t.PonS(1:numPanels, 2:numPanels) = C_tComb;

%% Generating systems of equations (3/4) - Primary
% The C matrices
[C_n1.P, C_n2.P, C_t1.P, C_t2.P] = deal(zeros(numPanels, numPanels));
                                     
% Control points on the primary, vortices also the primary
xC = xC_P;
xB = xB_P;
yC = yC_P;
yB = yB_P;
cosTheta = cosThetaP;
sinTheta = sinThetaP;
S_panel = S_panelP;
theta = thetaP;

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
        C_n2.P(i, j) = D + 0.5*(Q*F)/S_panel(j) - ...
                         (A*C + D*E)*G/S_panel(j);
        C_n1.P(i, j) = 0.5*D*F + C*G - C_n2.P(i, j);
        
        C_t2.P(i, j) = C + 0.5*(P*F)/S_panel(j) + ...
                         (A*D - C*E)*G/S_panel(j);
        C_t1.P(i, j) = 0.5*C*F - D*G - C_t2.P(i, j);
    end
end

% Fixing up the C matrices for cases where i = j (diagonal)
repDiagC_n1 = diag(ones(1, numPanels))*(-1); % replacement diagonal matrix
repDiagC_n2 = diag(ones(1, numPanels));
% Replacing diagonal of C_n1
C_n1.P = C_n1.P - diag(diag(C_n1.P)) + repDiagC_n1; 
C_n2.P = C_n2.P - diag(diag(C_n2.P)) + repDiagC_n2;

% Replacing diagonal of C_t1
repDiagC_t = diag(ones(1, numPanels))*(pi/2); % same replacement for 1 & 2
C_t1.P = C_t1.P - diag(diag(C_t1.P)) + repDiagC_t; 
C_t2.P = C_t2.P - diag(diag(C_t2.P)) + repDiagC_t;

% ******* Condensing the C matrices into A matrices. *******
% A_n matrix
A_n.P = zeros(numBounds, numBounds);
A_n.P(numBounds, 1) = 1;
A_n.P(numBounds, numBounds) = 1;
A_n.P(1:numPanels, 1) = C_n1.P(:, 1);
A_n.P(1:numPanels, numBounds) = C_n2.P(:, numPanels);
% Combining innards of C matrices as specified in Kuethe and Chow
C_nComb = C_n1.P(:, 2:numPanels) + C_n2.P(:, 1:numPanels - 1);
A_n.P(1:numPanels, 2:numPanels) = C_nComb; % voila

% A_t matrix
A_t.P = zeros(numPanels, numBounds);
A_t.P(:, 1) = C_t1.P(:, 1);
A_t.P(:, numBounds) = C_t2.P(:, numPanels);
% Combinign C_t matrices
C_tComb = C_t1.P(:, 2:numPanels) + C_t2.P(:, 1:numPanels - 1);
A_t.P(1:numPanels, 2:numPanels) = C_tComb;

%% Generating systems of equations (4/4) - Secondary
% The C matrices
[C_n1.S, C_n2.S, C_t1.S, C_t2.S] = ...
                                         deal(zeros(numPanels, numPanels));
                                     
% Control points on the secondary, vortices also the secondary
xC = xC_S;
xB = xB_S;
yC = yC_S;
yB = yB_S;
cosTheta = cosThetaS;
sinTheta = sinThetaS;
S_panel = S_panelS;
theta = thetaS;

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
        C_n2.S(i, j) = D + 0.5*(Q*F)/S_panel(j) - ...
                         (A*C + D*E)*G/S_panel(j);
        C_n1.S(i, j) = 0.5*D*F + C*G - C_n2.S(i, j);
        
        C_t2.S(i, j) = C + 0.5*(P*F)/S_panel(j) + ...
                         (A*D - C*E)*G/S_panel(j);
        C_t1.S(i, j) = 0.5*C*F - D*G - C_t2.S(i, j);
    end
end

% Fixing up the C matrices for cases where i = j (diagonal)
repDiagC_n1 = diag(ones(1, numPanels))*(-1); % replacement diagonal matrix
repDiagC_n2 = diag(ones(1, numPanels));
% Replacing diagonal of C_n1
C_n1.S = C_n1.S - diag(diag(C_n1.S)) + repDiagC_n1; 
C_n2.S = C_n2.S - diag(diag(C_n2.S)) + repDiagC_n2;

% Replacing diagonal of C_t1
repDiagC_t = diag(ones(1, numPanels))*(pi/2); % same replacement for 1 & 2
C_t1.S = C_t1.S - diag(diag(C_t1.S)) + repDiagC_t; 
C_t2.S = C_t2.S - diag(diag(C_t2.S)) + repDiagC_t;

% ******* Condensing the C matrices into A matrices. *******
% A_n matrix
A_n.S = zeros(numBounds, numBounds);
A_n.S(numBounds, 1) = 1;
A_n.S(numBounds, numBounds) = 1;
A_n.S(1:numPanels, 1) = C_n1.S(:, 1);
A_n.S(1:numPanels, numBounds) = C_n2.S(:, numPanels);
% Combining innards of C matrices as specified in Kuethe and Chow
C_nComb = C_n1.S(:, 2:numPanels) + C_n2.S(:, 1:numPanels - 1);
A_n.S(1:numPanels, 2:numPanels) = C_nComb; % voila

% A_t matrix - Also called the influence coefficient matrix. 
A_t.S = zeros(numPanels, numBounds);
A_t.S(:, 1) = C_t1.S(:, 1);
A_t.S(:, numBounds) = C_t2.S(:, numPanels);
% Combinign C_t matrices
C_tComb = C_t1.S(:, 2:numPanels) + C_t2.S(:, 1:numPanels - 1);
A_t.S(1:numPanels, 2:numPanels) = C_tComb;

%% Concatenating all A_n matrices to solve big boy system of eqns. 
A_nTot = zeros(2*numBounds, 2*numBounds);
A_nTot(1:numBounds, 1:numBounds) = A_n.P;
A_nTot(1:numBounds, numBounds + 1: 2*numBounds) = A_n.SonP;
A_nTot(numBounds+1:2*numBounds, 1:numBounds) = A_n.PonS;
A_nTot(numBounds+1:2*numBounds, numBounds+1:2*numBounds) = A_n.P;

RHStot = [RHS_P; RHS_S];
gama = A_nTot\RHStot;

gamma.P = gama(1:numBounds);
gamma.S = gama(numBounds + 1: 2*numBounds);
%% Determining Velocity using gamma

Vel_cont.P = zeros(1, numPanels);
Vel_cont.S = zeros(1, numPanels);
for i = 1:numPanels
    Vel_cont.P(i) = cos(thetaP(i) - alpha);
    Vel_cont.S(i) = cos(thetaS(i) - alpha);
    for j = 1:numBounds
        % Includes free stream, influences of the body's own panel
        % vortices, and the influence of the other body's panel vortices. 
        Vel_cont.P(i) = Vel_cont.P(i) + A_t.P(i, j)*gamma.P(j) + ...
                        A_t.SonP(i, j)*gamma.S(j); 
        Vel_cont.S(i) = Vel_cont.S(i) + A_t.S(i, j)*gamma.S(j) + ...
                        A_t.PonS(i, j)*gamma.P(j); 
    end
end

% Coefficient of pressure (Vel_cont is dimensionless)
c_pVec.P = 1 - (Vel_cont.P).^2;
c_pVec.S = 1 - (Vel_cont.S).^2;
%% Calculating circulation from velocity and panel length
GamTotP = sum(Vel_cont.P.*S_panelP);
GamTotS = sum(Vel_cont.S.*S_panelS);
cPrimary = xB_P(1);
cSecondary = (xB_S(1)^2 + yB_S(1)^2)^(1/2);
varargout{1} = 2*GamTotP/cPrimary + 2*GamTotS/cSecondary;


c_pVecTot = [c_pVec.P, c_pVec.S];
xCTot = [xC_P, xC_S];

if nargout == 2
    varargout{2} = [xCTot', c_pVecTot'];
end

end












