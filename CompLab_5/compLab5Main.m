%--------------------------------------------------------------------------
% compLab5Main perfoms the calculations required for Problems 2 and 3 for
% \textit{Computational Lab 5: Flow Over Finite Wings} (fancy italics for 
% fancy lab name). These inlude error analysis in relation to number of 
% Fourier terms used in Prandtl solution, as well as an analysis of the 
% relationship between span efficiency factor and aspect ratio.
%
% Inputs:
%        None
% Outputs: 
%       - A plot detailing error trends with increased Fourier terms in 
%         PLLT solution (log and linear scales).
%       - A plot detailing relationship between span efficiency factor and
%         aspect ratio.
% Dependancies:
%       - prandtlLLT.m
%       - NACAairfoilPlot.m
%       - vortexPanel.m
%       - hline.m
%
% Created: 11/3/17 - Connor Ott
% Last Modified: 11/9/17 - Connor Ott
%--------------------------------------------------------------------------

clear; close all; clc;

%% Getting lift slope for root and tip airfoil shapes.
% Determining error using a wing with NACA 4412 at root and 2412 at tip.
% Need lift slope -> Using NACAairfoilPlot and vortexPanel

% The next 50 or so lines are directly from my own comp lab 4 code, this is
% me citing myself.

NACA = ['4412'; '2412'];

alphaMax = -5;  % [deg]
alphaMin = 5;   % [deg]
numAlpha = 15;
alphaVec = linspace(alphaMin, alphaMax, numAlpha);
[numAfoil, ~] = size(NACA);

M = str2num(NACA(:, 1));
P = str2num(NACA(:, 2));
t = str2num(NACA(:, 3:4));
c = 1;      % [m] chord length
V_inf = 0;  % [m/s] - This variable has no effect

[xPts, yPts, cpPlot] = deal(cell(1, numAfoil));
clData = zeros(numAlpha, numAfoil);
N = 50; % Number of panels used for determining cl at each alpha
for i = 1:numAfoil
    for j = 1:numAlpha
        % Getting C_l data
        [xPts{i}, yPts{i}] = NACAairfoilPlot(M(i), P(i), t(i), c, N);
        clData(j, i) = vortexPanel(xPts{i}, yPts{i}, V_inf, alphaVec(j));
    end
end

% Lift slopes for the root and chord from polyfit of alpha-cl plots.
alpha_cl{1} = polyfit(alphaVec', clData(:, 1), 1);
alpha_cl{2} = polyfit(alphaVec', clData(:, 2), 1);

a0_r = alpha_cl{1}(1) * 180/pi; % [rad^-1]
a0_t = alpha_cl{2}(1) * 180/pi;

% Swapping cl and alpha data to get 'alpha' intercept from polyfit. 
% Translates to zero-lift angle of attack.
cl_alpha{1} = polyfit(clData(:, 1), alphaVec', 1);
cl_alpha{2} = polyfit(clData(:, 2), alphaVec', 1);

aero_r = cl_alpha{1}(2); % [deg]
aero_t = cl_alpha{2}(2);

% From Problem 2 Statement:
geo_r   = 7;  % [deg]
geo_t   = 4;  % [deg]
c_r     = 6;  % [ft]
c_t     = 2;  % [ft]
b       = 40; % [ft]

% Calculating error for wing defined above.
N_exact = 300; % Terms in circulation solution from PLTT
[e_exact, C_LExact, C_DiExact] = prandtlLLT(b, a0_t, a0_r, c_t, c_r, ...
                                            aero_t, aero_r, geo_t, ...
                                            geo_r, N_exact);

MPH2FTPS = 5280/3600;   % Conversion factor                                   
V_inf = 130 * MPH2FTPS; % [ft/s]
rho_inf = 2.3769e-3;    % [slug/ft^3] - standard atmosphere table
q_inf = 0.5 * rho_inf * V_inf.^2;  % [Pa]
S = b/2 * (c_t + c_r);             % [ft^2]
L_exact = C_LExact * q_inf * S;    % [lbs]
Di_exact = C_DiExact * q_inf * S;  % [lbs]

trialVec = 3:40;
errMat = zeros(2, length(trialVec));
j = 1;
k = 1;
for i = 1:length(trialVec)
    N_trial = trialVec(i);
    [e_trial, C_LTrial, C_DiTrial] = prandtlLLT(b, a0_t, a0_r, c_t, ...
        c_r, aero_t, aero_r, geo_t, geo_r, N_trial);
    
    L_trial = C_LTrial * q_inf * S;    % [lbs]
    Di_trial = C_DiTrial * q_inf * S;  % [lbs]
    
    errMat(1, i) = abs(L_trial - L_exact) / L_exact * 100;
    errMat(2, i) = abs(Di_trial - Di_exact) / Di_exact * 100;
end

nomErrVec = [5, 1, 0.1];
for i = 1:3
    logicTempL = errMat(1, :) <= nomErrVec(i);
    tempidxL = find(logicTempL);
    tempNL = trialVec(tempidxL(1));
    nomLift(i) = tempNL;
    
    logicTempD = errMat(2, :) <= nomErrVec(i);
    tempidxD = find(logicTempD);
    tempND = trialVec(tempidxD(1));
    nomDrag(i) = tempND;
end

T = table(nomErrVec', nomDrag', nomLift', ...
    'VariableNames', ...
    {'Error_Percent',...
    'Drag_Term_Num', ...
    'Lift_Term_Num'});
disp(T)
fprintf('Exact Lift Value: <strong>%.1f lbs </strong> \n', L_exact);
fprintf('Exact Induced Drag Value: <strong>%.1f lbs </strong> \n\n', ...
         Di_exact);
%% Plotting and fitting error
% Definitely doing an error fit again
errFit.c_L = fit(trialVec', errMat(1, :)', 'exp2');
errFit.c_Di = fit(trialVec', errMat(2, :)', 'exp2');
fitNVec = linspace(1, 50, 1000);
errFitVec(1, :) = errFit.c_L(fitNVec);
errFitVec(2, :) = errFit.c_Di(fitNVec);

figure
set(0, 'defaulttextinterpreter', 'latex')

% Drag Plot
subplot(1, 2, 1)
hold on; grid on; grid minor;

plot(trialVec, errMat(2, :), 'ob', 'linewidth', 1)
plot(fitNVec, errFitVec(2, :), 'b-.', 'linewidth', 0.75)

hline(5, 'k--', '$5 \%$ Error');
hline(1, 'k--', '$1 \%$ Error');
hline(0.1, 'k--', '$0.1\%$ Error');

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 13, ...
         'box', 'on'); 
axis([1 trialVec(end)*1.2 0 max(max(errMat))*1.2])
xlabel('Number of Terms Used')   
ylabel(sprintf('Percent Errror Relative to %.f Terms', N_exact))
leg = legend('Induced Drag Error', ...
              sprintf(['Drag Error Fit: E = $%.1fe^{%.1fN}$ + ', ...
                       '$%.1fe^{%.1fN}$'], ...
                       errFit.c_Di.a, errFit.c_Di.b, errFit.c_Di.c, ...
                       errFit.c_Di.d));
set(leg, 'Interpreter', 'latex',...
         'fontsize', 11);
     
% Lift Plot
subplot(1, 2, 2)
hold on; grid on; grid minor;

plot(trialVec, errMat(1, :), 'sr', 'linewidth', 1) 
plot(fitNVec, errFitVec(1, :), 'r-.', 'linewidth', 0.75)

hline(5, 'k--', '$5 \%$ Error');
hline(1, 'k--', '$1 \%$ Error');
hline(0.1, 'k--', '$0.1\%$ Error');

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 13, ...
         'box', 'on'); 
axis([1 trialVec(end)*1.2 0 max(max(errMat))*1.2])
xlabel('Number of Terms Used')   
ylabel(sprintf('Percent Errror Relative to %.f Terms', N_exact))
leg = legend('Lift Error', ... 
             sprintf(['Lift Error Fit: E = $%.1fe^{%.1fN}$ + ', ...
                      '$%.1fe^{%.1fN}$'], ...
                      errFit.c_L.a, errFit.c_L.b, errFit.c_L.c, ...
                      errFit.c_L.d));
set(leg, 'Interpreter', 'latex',...
         'fontsize', 12);
     
% Super title - goes over both subplots
ha = axes('Position',[0 0 1 1], ...
          'Xlim',[0 1], 'Ylim',[0 1], ...
          'Box', 'off', ...
          'Visible', 'off', ...
          'Units', 'normalized',  ...
          'clipping', 'off');
text(0.5, 0.99,'Relative Error vs. Number of Terms',...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'top', ...
               'fontsize', 17);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

saveas(gcf, 'LDerror.png')

%% Log-Log plot
figure
% Lift Plot
subplot(1, 2, 1)
hold on; grid on; grid minor;

plot(trialVec, errMat(2, :), 'ob-', 'linewidth', 1)

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 13, ...
         'box', 'on', ...
         'xscale', 'log', 'yscale', 'log'); 
axis([2 trialVec(end)*1.2 0 max(errMat(2, :))*1.2])
xlabel('Number of Terms Used')   
ylabel(sprintf('Percent Errror Relative to %.f Terms', N_exact))
leg = legend('Induced Drag Error');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 11);
     
% Drag Plot
subplot(1, 2, 2)
hold on; grid on; grid minor;

plot(trialVec, errMat(1, :), 'sr-', 'linewidth', 1) 

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 13, ...
         'box', 'on', ...
         'xscale', 'log', 'yscale', 'log');
axis([2 trialVec(end)*1.2 0 max(errMat(1, :))*1.2])
xlabel('Number of Terms Used')   
ylabel(sprintf('Percent Errror Relative to %.f Terms', N_exact))
leg = legend('Lift Error');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 12);
     
% Super title - Goes over both subplots
ha = axes('Position',[0 0 1 1], ...
          'Xlim',[0 1], 'Ylim',[0 1], ...
          'Box', 'off', ...
          'Visible', 'off', ...
          'Units', 'normalized',  ...
          'clipping', 'off');
text(0.5, 0.99,'Relative Error vs. Number of Terms',...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'top', ...
               'fontsize', 17);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

saveas(gcf, 'LDerrorLOG.png')


%% Aspect Ratio and c_t/c_r vs. Span Efficiency
a0_t = 2*pi;
a0_r = 2*pi;
aero_t = 7;
aero_r = 7;
geo_t = 0;
geo_r = 0;
N = 20;
S = 240; % [ft^2] 

% Keeping S constant and varying AR, taper ratio, c_t, and c_r.
lamVec = linspace(0, 1, 250); % c_t/c_r ratio
AR_vec = [4, 6, 8, 10];
bVec = sqrt(S .* AR_vec);
% not to be confused with errMat, eMat holds span efficiency values.
eMat = zeros(length(bVec), length(lamVec)); 

for i = 1:length(AR_vec)
    for j = 1:length(lamVec)
        c_r =  (2 * S) / (bVec(i) * (lamVec(j) + 1));
        c_t = c_r * lamVec(j);
        b = bVec(i);
        [eMat(i, j), ~] = prandtlLLT(b, a0_t, a0_r, c_t, c_r, aero_t, ...
                                     aero_r, geo_t, geo_r, N);
    end
    eMax(i) = max(eMat(i, :));
    logLam = eMat(i, :) == eMax(i);
    maxELoc = find(logLam);
    lamEmax(i) = lamVec(maxELoc);
end
optLamMean = mean(lamEmax);
T = table(lamEmax', eMax', AR_vec', ...
    'variablenames', ...
    {'Taper_Ratio', ...
    'Span_efficiency', ...
    'Aspect_Ratio'});
disp(T)
fprintf('Mean Optimal Taper Ratio: <strong>%.4f </strong> \n', optLamMean)

%% Plotting
figure
hold on; grid on; grid minor;
axis([0, 1, 0.84, 1])

plot(lamVec, eMat(1, :), 'k-', 'linewidth', 1.2)
plot(lamVec, eMat(2, :), 'b-', 'linewidth', 1.2)
plot(lamVec, eMat(3, :), 'm-', 'linewidth', 1.2)
plot(lamVec, eMat(4, :), 'c-', 'linewidth', 1.2)

for i = 1:length(AR_vec)
    plot(lamEmax(i), eMax(i), 'rx', 'linewidth', 1) 
end

set(gca, 'TickLabelInterpreter', 'latex',...
         'fontsize', 11, ...
         'box', 'on')
title(sprintf(['Taper Ratio vs. Span Efficiency for ', ...
              'Various Aspect Ratios']))
xlabel('Taper Ratio, $\lambda$')   
ylabel('Span Efficiency Factor, $e$')
leg = legend(sprintf('AR = %.f', AR_vec(1)), ... 
             sprintf('AR = %.f', AR_vec(2)), ...
             sprintf('AR = %.f', AR_vec(3)), ...
             sprintf('AR = %.f', AR_vec(4)), ...
             'Maximum e For Given AR', ...
             'location', 'best');
set(leg, 'Interpreter', 'latex',...
         'fontsize', 9);
saveas(gcf, 'AR_taper_e.png')

