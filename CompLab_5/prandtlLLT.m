%--------------------------------------------------------------------------
% prandtlLLT predicts span efficiency factor, lift coefficient, and induced
% drag coefficient for a finite wing. 
%
% [e, C_L, C_Di] = prandtlLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r,...
%                             geo_t, geo_r, N)
%
% Inputs: (Use consistent units)
%       b         - Span 
%       a0        - [rad^-1] Lift slope
%       c         - Chord length
%       aero      - [deg] zero-lift angle of attack
%       geo       - [deg] geometric angle of attack
%
%       _t        - Denotes wing tip characteristic
%       _r        - Denotes wing root characteristic
%
% Outputs:
%       e         - Span efficiency factor
%       C_L       - Lift Coefficient
%       C_Di      - Induced Drag Coefficient
%
% Created: 11/3/17 - Connor Ott
% Last Modified: 11/7/17 - Connor Ott
%--------------------------------------------------------------------------
function [e, C_L, C_Di] = prandtlLLT(b, a0_t, a0_r, c_t, c_r, aero_t, ...
                                     aero_r, geo_t, geo_r, N)
                            
% Taking the form h = m(theta) + b where h is an airfoil property.
m_a0 = (a0_t - a0_r);
m_c = (c_t - c_r);
m_aero = (aero_t - aero_r);
m_geo = (geo_t - geo_r);

% Defining theta = 0 at root and theta = pi/2 at tip.
aSlope = @(thta) m_a0*cos(thta) + a0_r;
AoA = @(thta) m_geo*cos(thta) + geo_r;
AoA_L0 = @(thta) m_aero*cos(thta) + aero_r;
chord = @(thta) m_c*cos(thta) + c_r;

AR = b^2 / ((c_r + c_t)/2 * b); % Used in coefficient calculation
 
A = zeros(N, N);
RHS = zeros(N, 1) ;
thetaInit = linspace(0, pi/2, N+1);
thetaTrim = thetaInit(2:end); % trimming an end points because it's bad.
for i = 1:N
    for j = 1:N
        n = 2*j - 1;
        theta = thetaTrim(i);
        % Prandtl LLT 
        A(i, j) = 4*b/(aSlope(theta)*(chord(theta)))*sin(n*theta) +...
                      (n*sin(n*theta)/(sin(theta)));
        
        RHS(i) = (AoA(theta) - AoA_L0(theta)) * pi/180;
    end
end

% Solving for Fourier coefficients
coefVec = A\RHS;

% Summing for delta to find span efficieny factor, e = (1 + delta)^(-1)
delta = 0;
for i = 2:N
    n = 2*i - 1;
    delta = delta + n*(coefVec(i)/coefVec(1))^2;
end

e = (delta + 1) ^-1;
C_L = coefVec(1) * pi * AR;
C_Di = C_L^2 / (pi*e*AR);      
end