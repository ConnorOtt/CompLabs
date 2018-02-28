% -------------------------------------------------------------------------
% 
% cl_vs_alphaAnalysis Plots sectional coefficients of lift vs. angle of 
% attack for various NACA series airoils 
%
%   Inputs: 
%       - None
%   Outputs:
%       - A cl-alpha diagram with data for each airfoil.
%       - A cp-x diagram for the NACA 0012.
%       - A plot detailing error analysis for numbers of panels. 
%   Dependancies: 
%       - NACAairfoilPlot.m
%       - vortexPanel.m
%       - hline.m
%       - vline.m
%       - actualVars.mat (locks inputs for error analysis)
%       - actualMats.mat (High N vortex panel results, for comparison in 
%                        error calculations)
%
% Created: 10/10/17 - Connor Ott
% Last Modified: 10/21/17 - Connor Ott
% -------------------------------------------------------------------------

clear; close all; clc;

NACA = ['0012'; '2412'; '4412'; '2430'];

alphaMax = 15;  % [deg]
alphaMin = -5; % [deg]
numAlpha = 15;
alphaVec = linspace(alphaMin, alphaMax, numAlpha);
[numAfoil, ~] = size(NACA);

M = str2num(NACA(:, 1));
P = str2num(NACA(:, 2));
t = str2num(NACA(:, 3:4));
c = 2;        % [m] chord length - reassigned throughout code frequently
V_inf = 50;   % [m/s] Free stream velocity

[xPts, yPts, cpPlot] = deal(cell(1, numAfoil));
clData = zeros(numAlpha, numAfoil);

%% Getting c_l data for each airfoil at each angle of attack 
N = 50; % Number of panels used for determining cl at each alpha
for i = 1:numAfoil
    for j = 1:numAlpha
        [xPts{i}, yPts{i}] = NACAairfoilPlot(M(i), P(i), t(i), c, N, ...
                                            'HalfCos');
        clData(j, i) = vortexPanel(xPts{i}, yPts{i}, V_inf, alphaVec(j));
    end
end

%% Example x vs. cp diagram
c = 2;
N = 500;
alpha = 8;
delta = 10;
[exX, exY] = NACAairfoilPlot(M(1), P(1), t(1), c, N, 'HalfCos');
[~, exCp] = vortexPanel(exX, exY, V_inf, alpha); 

figure
set(0, 'defaulttextinterpreter', 'latex')
hold on; grid on; grid minor;
axis([0 c min(exCp(:,2))*1.1 max(exCp(:,2))*1.1])

% Emphasizing c_p = 0
CpZero = hline(0, 'k','$c_p$ = 0');
set(CpZero, 'handlevisibility','off', ...
         'color', [0.85, 0.2, 0.2], ...
         'lineWidth', 1.1);

plot(exCp(:, 1), exCp(:, 2), 'xb--');
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
title('Sectional Pressure Coefficient vs. x Position Along Chord')
xlabel('x Position Along Chord [m]')   
ylabel('$C_p$')
set(gca, 'YDir', 'reverse')
saveas(gcf, 'cpPlot.png');


%% Plotting angle of attack vs. alpha for each airfoil
figure
hold on; grid on; grid minor;
axis([alphaMin*1.2, alphaMax*1.1, min(min(clData))*1.1, ...
      max(max(clData))*1.05])

% Identifying 0 alpha and 0 c_l
ClZero = hline(0, 'k','$c_l$ = 0');
yAx = vline(0, 'k', '$\alpha$ = $0^{\circ}$');
set(ClZero, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);
set(yAx, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);
     
plot(alphaVec, clData(:, 1) ,'--rs', 'linewidth' ,0.5);
plot(alphaVec, clData(:, 2) ,'--bo', 'linewidth' ,0.5);
plot(alphaVec, clData(:, 3) ,'--mv', 'linewidth' ,0.5);
plot(alphaVec, clData(:, 4) ,'--k+', 'linewidth' ,0.5);

set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
title('Sectional Lift Coefficient vs. Angle of Attack')
xlabel('Angle of Attack, $\alpha$ $[^{\circ}]$')   
ylabel('Sectional Lift Coefficient, $c_l$')
leg = legend([NACA(1, :), ' Lift Coefficient'], ...
             [NACA(2, :), ' Lift Coefficient'], ...
             [NACA(3, :), ' Lift Coefficient'], ...
             [NACA(4, :), ' Lift Coefficient'], ...
             'location', 'northwest');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 10);
saveas(gcf, 'clPlot.png');


%% Looking into dcl/dalphas
for i = 1:numAfoil
   tempfit = polyfit(alphaVec' ,clData(:, i), 1) * 180/pi;
   dcl_dalpha(i) = tempfit(1); 
   tempfit2 = polyfit(clData(:, i), alphaVec', 1);
   alphaL0(i) = tempfit2(2);

end

%% Looking into Error for different numbers of panels.
% Only using a NACA 0012

% Loading data and variables fixed to the accepted Cp curve, used to
% compare trials of different N values to the accepted
load('actualMats.mat')
load('actualVars.mat')

% Logorithmically spaced trials
minTrial = 5e2;       % Minimum panels used in error calc 
maxTrial = 4e3;       % Maximum panels used in error calc
numTrials = 8;        % Number of panel numbers used in error calc
% trialVec = exp(linspace(log(minTrial), log(maxTrial), numTrials));
trialVec = nonLinspace(minTrial, maxTrial, numTrials, 'exp10');
[eMatLin, eMatCos] = deal(zeros(1, numTrials));
%% Evaluating and comparing at different N values

for i = 1:numTrials

    nTrial = floor(trialVec(i));
    [xBLin, yBLin] = NACAairfoilPlot(0, 0, 12, c, nTrial);
    [xBCos, yBCos] = NACAairfoilPlot(0, 0, 12, c, nTrial, 'HalfCos');   
    
    [~, cPLin] = vortexPanel(xBLin, yBLin, V_inf, alpha);
    [~, cPCos] = vortexPanel(xBCos, yBCos, V_inf, alpha);
    
    % Interpolating between actual panel location
    % x points must be unique so upper and lower surfaces must be done
    % seperately.
    cpInterpLin.low = interp1(cPLin(1:end/2, 1), cPLin(1:end/2, 2), ...
                             cPActLin(1:end/2, 1), ...
                             'linear', ...
                             'extrap');
    cpInterpLin.up = interp1(cPLin(end/2+1:end, 1), ...
                             cPLin(end/2+1:end, 2), ...
                             cPActLin(end/2+1:end,1) ,...
                             'linear', ...
                             'extrap');
    cpInterpLin.tot = [cpInterpLin.low; cpInterpLin.up];
    
    cpInterpCos.low = interp1(cPCos(1:end/2, 1), cPCos(1:end/2, 2), ...
                             cPActCos(1:end/2, 1), ...
                             'linear', ...
                             'extrap');
    cpInterpCos.up = interp1(cPCos(end/2+1:end, 1), ...
                             cPCos(end/2+1:end, 2), ...
                             cPActCos(end/2+1:end,1) ,...
                             'linear', ...
                             'extrap');
    cpInterpCos.tot = [cpInterpCos.low; cpInterpCos.up];  
    
    % Percent error defined as 2-norm of difference between the trial 
    % vector and accepted vector / 2-norm of accepted vector * 100. 
    eMatLin(i) = norm(cpInterpLin.tot - cPActLin(:, 2))/...
                 norm(cPActLin(:, 2))*100;
    eMatCos(i) = norm(cpInterpCos.tot - cPActCos(:, 2))/...
                 norm(cPActCos(:, 2))*100;
    trialVec(i) = nTrial;
end


%% Error plotting

% Solid HELL YEAH on the plot fitting
errorFitLin = fit(trialVec', eMatLin', 'power2');
errorFitCos = fit(trialVec', eMatCos', 'power2');
fitVec = linspace(1, trialVec(end)*1.2, 1e4);
errorFitVecCos = errorFitCos(fitVec); 
errorFitVecLin = errorFitLin(fitVec);

figure
hold on; grid on; grid minor;
plot(trialVec, eMatLin, 'sb', 'linewidth', 0.75)
plot(trialVec, eMatCos, 'or', 'linewidth', 0.75) 
plot(fitVec, errorFitVecLin, 'b-.', 'linewidth', 1)
plot(fitVec, errorFitVecCos, 'r-.', 'linewidth', 1)

axis([0 trialVec(end)*1.2 0 max(eMatLin)*1.1])

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 11)
title(sprintf('Relative Error vs. Number of Panels'))
xlabel('Number of Panels Used')   
ylabel('$\%$ Error Relative to $10^{4}$ Panels')
leg = legend('Linear Spaced Panel Error', ... 
             'Cosine-Spaced Panel Error', ...
             sprintf(['Error Trend Fit (Linear) $E = %.2f N^{%.2f}', ...
             '%.2f$'], errorFitLin.a, errorFitLin.b, errorFitLin.c), ...
             sprintf(['Error Trend Fit (Cosine) $E = %.2f N^{%.2f}', ...
             '+%.2f$'], errorFitCos.a, errorFitCos.b, errorFitCos.c));
set(leg, 'Interpreter', 'latex',...
         'fontsize', 9);
saveas(gcf, 'errorTrend.png');
     

%% Nominal number of panels
Nnom = 969; 

alphaVecNom = [-5, 0, 5, 10];
[xPtsNom, yPtsNom] = NACAairfoilPlot(M(1), P(1), t(1), c, Nnom, ...
                                     'HalfCos');
for i = 1:length(alphaVecNom)
    for j = 1:length(alphaVecNom)
        [clData2(j, i), cpData{j}] = vortexPanel(xPtsNom, yPtsNom,...
                                     V_inf, alphaVecNom(j));
    end
end

%% Plotting and such
figure
hold on; grid on; grid minor;
set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 11)
title('Sectional Pressure Coefficient vs. x Position Along Chord')
xlabel('x Position Along Chord [m]')   
ylabel('$C_p$')
set(gca, 'YDir', 'reverse')

plot(cpData{1}(:, 1), cpData{1}(:, 2) ,'r', 'linewidth', 2);
plot(cpData{2}(:, 1), cpData{2}(:, 2) ,'b', 'linewidth', 2);
plot(cpData{3}(:, 1), cpData{3}(:, 2) ,'c', 'linewidth', 2);
plot(cpData{4}(:, 1), cpData{4}(:, 2) ,'k', 'linewidth', 2);

leg = legend(['Pressure Coefficient for $\alpha = ', ...
              num2str(alphaVecNom(1)), '^{\circ}$'], ...
             ['Pressure Coefficient for $\alpha = ', ...
              num2str(alphaVecNom(2)), '^{\circ}$'], ...
             ['Pressure Coefficient for $\alpha = ', ...
              num2str(alphaVecNom(3)), '^{\circ}$'], ...
             ['Pressure Coefficient for $\alpha = ', ...
              num2str(alphaVecNom(4)), '^{\circ}$'], ...
             'location', 'northeast');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 9);
saveas(gcf, 'alphaNomComp.png')

figure;
hold on; grid on; grid minor;
axis([alphaVecNom(1)*1.2, alphaVecNom(end)*1.1, min(min(clData2))*1.1, ...
      max(max(clData2))*1.05])
% Identifying 0 alpha and 0 c_l
ClZero = hline(0, 'k','$c_l$ = 0');
yAx = vline(0, 'k', '$\alpha$ = $0^{\circ}$');
set(ClZero, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);
set(yAx, 'handlevisibility','off', ...
         'color', [0, 0, 0], ...
         'lineWidth', 0.5);

plot(alphaVecNom, clData2(:, 1) ,'--rs', 'linewidth' ,0.5);
     
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
title('Sectional Lift Coefficient vs. Angle of Attack')
xlabel('Angle of Attack, $\alpha$ $[^{\circ}]$')   
ylabel('Sectional Lift Coefficient, $c_l$')
leg = legend([NACA(1, :), ' Lift Coefficient'], ...
             'location', 'northwest');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 10);
saveas(gcf, 'alphaNomComp_cl.png')
