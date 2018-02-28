% *************************************************************************
%
% flapVnoflap Compares sectional coefficient of lift of a NACA 0012 airfoil
% with and without a trailing edge flap. 
%
% Dependancies:
%      - vortexPanel.m
%      - vortexPanel_flap.m
%      - hline.m (plotting)
%      - vline.m (plotting)
%
% Created: 10/20/17 - Connor Ott
% Last Modified: 10/10/17 - Connor Ott
%
% *************************************************************************
clc; clear; close all

M = 0;
P = 0;
t = 12;
c = 2;          % [m] Chord length of airfoil
N = 40;         % Number of panels
delta = 10;     % [deg] Deflection angle of flap 
V_inf = 100;    % [m/s] Fee stream velocity ()
alpha = 5;      % [deg] Angle of attack of airfoil

% getting positions for flap and airfoils
[x1, y1] = NACAairfoilPlot(M, P, t, c, N, 'HalfCos');
[x2, y2] = NACAairfoilPlot(M, P, t, c*0.35, N, ...
                           'HalfCos', ...
                           'delta', delta);

% Shifting flap back behind the airfoil
x2 = x2 + c + 0.05*c;

% plot(x1, y1, x2, y2)
% axis equal   
[clFlap, cPFlap] = vortexPanel_flap(x1, y1, x2, y2, V_inf, alpha);
[clNoFlap, cPNoFlap] = vortexPanel(x1, y1, V_inf, alpha);

%% Cp comparison
figure
set(0, 'defaulttextinterpreter', 'latex')
hold on; grid on; grid minor;
axis([-cPFlap(end,1)*.05,  cPFlap(end,1)*1.05,... 
      min(cPFlap(:, 2))*1.1, max(cPFlap(:, 2))*1.1])

% Emphasizing c_p = 0
CpZero = hline(0, 'k','$c_p$ = 0');
set(CpZero, 'handlevisibility','off', ...
         'color', [0.85, 0.2, 0.2], ...
         'lineWidth', 1.1);

plot(cPNoFlap(:, 1), cPNoFlap(:, 2), 'xr--', 'linewidth', 0.5);
plot(cPFlap(1:N, 1), cPFlap(1:N, 2), 'xb--', 'linewidth', 0.5);
plot(cPFlap(N+1:2*N, 1), cPFlap(N+1:2*N, 2), 'xb--', 'linewidth', 0.5);

set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
title('Sectional Pressure Coefficient vs. x Position Along Chord')
xlabel('x Position Along Chord [m]')   
ylabel('$C_p$')
set(gca, 'YDir', 'reverse')
leg = legend('$C_p$ Curve w/o Flap', ... 
             '$C_p$ Curve with Flap');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 9);
saveas(gcf, 'cpCompFlap.png');


%% Lift Slop Comparison
numAlpha = 15;
alphaMin = -10;
alphaMax = 20;
alphaVec = linspace(alphaMin, alphaMax, numAlpha);

clData = zeros(numAlpha, 2);

for i = 1:numAlpha
    clData(i, 1) = vortexPanel(x1, y1, V_inf, alphaVec(i));
    clData(i, 2) = vortexPanel_flap(x1, y1, x2, y2, V_inf, alphaVec(i));
end

figure
hold on; grid on; grid minor;
axis([min(alphaVec) max(alphaVec) min(min(clData)) max(max(clData))]);

CpZero = hline(0, 'k','$c_l$ = 0');
yAx = vline(0, 'k', '$\alpha$ = $0^{\circ}$');
set(CpZero, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);
set(yAx, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);
     
plot(alphaVec, clData(:, 1), 'rs--');
plot(alphaVec, clData(:, 2), 'bo--');

set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
title(sprintf('Angle of Attack Vs. $C_l$'))
xlabel('Angle of Attack, $\alpha~[^{\circ}]$')   
ylabel('Sectional Lift Coefficient, $C_l$')
leg = legend('Lift Slope w/o Flap','Lift Slope w/ Flap');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 10, ...
         'Location', 'NorthWest');
saveas(gcf, 'clCompFlap.png');







