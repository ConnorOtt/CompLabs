% -------------------------------------------------------------------------
% [Xpts, Ypts] = airfoilPlot(M, P, t, c, N) gives x and y coordinates 
% corresponding to the surface of any 4-digit NACA series airfoil.
% 
% [Xpts, Ypts] = airfoilPlot(M, P, t, c, N, 'delta', delta) Specifies a
% deflection angle, delta, in degrees.
%
% [Xpts, Ypts] = airfoilPlot(M, P, t, c, N, 'HalfCos') Specifies half
% cosing spacing, putting more panels towards the leading and trailing 
% egdes of the airfoil (default is equal-length panels). Credit to Divahar 
% Jayaraman (j.divahar@yahoo.com) for the half-cosine spacing functionality
% found in naca4gen.m (included in this folder).
%
% Inputs:
%       M -     Max Camber Value [% of chord]
%       P -     Max Camber Location {*10% of chord}
%       t -     Max thickness [% of chord]
%       c -     Chord length [m]
%       N -     Number of panels requested for vortex sheet simulation. Use
%               an even number so that there are an equal number of panels 
%               above and below the airfoil. If an odd N is input, N + 1
%               panels will be used for simplicity
%       delta - Angle of attack of airfoil. 
%
% Outputs:
%       xPts -  x coordinates of surface of airfoil, beginning and ending
%               at the trailing edge.
%       yPts -  y coordinates of the surface of the airfoil, also beginning
%               and ending at trailing edge. 
%
% Created: 10/06/17 - Connor Ott
% Last Modified: 10/19/17 - Connor Ott
% -------------------------------------------------------------------------


function [xPts, yPts] = NACAairfoilPlot(M, P, t, c, N, varargin)

delta = 0;
halfCos = 0;
if nargin > 5
    for i = 1:nargin - 5
        switch varargin{i}
            case 'delta'
                delta = varargin{i + 1};
            case 'HalfCos'
                halfCos = 1;
        end
    end
end

M = M/100;
P = P/10;
t = t/100;
if mod(N, 2) ~= 0 % Should be an even number of panels
   N = N + 1; % eat me 
end
NptsTop = N/2 + 1; % number of boundaries defining panels

[xPts, yPts] = deal(zeros(1, NptsTop));
%% Creating X and Y coordinates for other x and y coordinates or something.
if halfCos
    % Cosine spacing (D. Jayaraman)
    xVecTotal = nonLinspace(0, c, NptsTop, 'cos'); 
else
    xVecTotal = linspace(0, c, NptsTop);
end

xMC_idx = find(xVecTotal <= c*P);  % max camber indices.
xVecFrnt = xVecTotal(1:xMC_idx(end));
xVecBack = xVecTotal(xMC_idx(end)+1: end);

%% Airfoil geometry equations
yTness = @(x) t/0.2 .* c .* (0.2969*(x/c).^(0.5) - 0.1260*(x/c) - ...
    0.3516*(x/c).^2 + 0.2843*(x/c).^3 - 0.1036*(x/c).^4);

if P == 0 % If the airfoil is symmetric, call it quits
    
    xPts = [fliplr(xVecTotal(2:end)), xVecTotal];
    yPts = [-fliplr(yTness(xVecTotal(2:end))), yTness(xVecTotal)];
    
else % lets play hardball
    
    % Camber
    yCam_frnt = @(x) M .* x./(P.^2) .* (2*P - x/c); 
    yCam_back = @(x) M .* (c-x)./((1-P).^2) .* (1 + x/c - 2*P);

    dyCamF_dx = @(x) 2*M/(c*P^2).*(c*P - x);
    dyCamB_dx = @(x) 2*M/(c*(1-P)^2).*(c*P - x);
    
    % More Camber geometry
    xsi_frnt = @(x) atan(dyCamF_dx(x));
    xsi_back = @(x) atan(dyCamB_dx(x));
    
    % Evaluating at discrete locations
    yThick = yTness(xVecTotal);
    yCam = [yCam_frnt(xVecFrnt), yCam_back(xVecBack)];
    xsi = [xsi_frnt(xVecFrnt), xsi_back(xVecBack)];
    
    % Generating x and y coordinates for the surface of the airfoil.
    X_up = xVecTotal - yThick.*sin(xsi);
    Y_up = yCam + yThick.*cos(xsi);
    
    X_low = xVecTotal + yThick.*sin(xsi);
    Y_low = yCam - yThick.*cos(xsi);
    
    % Making data go from trailing edge around and back to trailing edge.
    xPts = [fliplr(X_low(2:end)), X_up];
    yPts = [fliplr(Y_low(2:end)), Y_up];
    
end

delta = delta * pi/180;
xPts = xPts*cos(-delta) - yPts*sin(-delta);
yPts = xPts*sin(-delta) + yPts*cos(-delta);
end