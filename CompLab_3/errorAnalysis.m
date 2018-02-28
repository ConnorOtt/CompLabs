% *************************************************************************
% 
% errorAnalysis - Determines N values required for different tolerances of
% accuracy. Also gives plots indicating areas of greatest uncertainty on
% the airfoil. 
%
% Inputs: 
%       -None
% Outputs:
%       -A plot showing varying error around the leading edge of the thin
%       airfoil.
%       -A plot showing varying error with increased number of vortices.
%       -Two print statements showing the fit lines for the data in the
%       plot second plot. 
% Dependancies:
%       -V_true.mat
%       -P_true.mat
%       -NOplotAirfoilFlow.m
%
% Created: 10/03/17 - Connor Ott
% Last Modified: 10/06/17 - Connor Ott
% 
% *************************************************************************


clc;clear;close all
tic

c = 2;              % [m] Chord Length
alpha = 10;         % [deg] Angle of Attack
V_inf = 100;        % [m/s] Free stream velocity
P_inf = 2.65e4;     % [Pa] Free stream pressure
rho_inf = 0.4135;   % [kg/m^3] Free stream density

%% Error Calculation
% Loading results from NOplotAirfoilFlow.m with 1e5 vortices (takes a long 
% time, don't want to spend time getting it on every run).
load('V_true.mat');
load('P_true.mat');

% Getting values at lower N for comparison
numErr = 7;
[P_N, V_N] = deal(cell(1, numErr));

% For holding error from each field.
[VerrVec, PerrVec, nVec] = deal(zeros(1, numErr)');

for n = 1:numErr
    % Num vortice pattern: 
    % 5e0 1e1 5e1 1e2 5e2 1e3 5e3 ... 
    if mod(n-1, 2) == 0 % Deciding number of vortices
        N = 5 * 10^((n-1)/2); 
    else
        N = 10^((n)/2); 
    end
    [P_N{n}.mat, V_N{n}.mat, P_N{n}.xMat, P_N{n}.yMat] = ...
        NOplotAirfoilFlow(c, alpha, V_inf, P_inf, rho_inf, N);
    
    % Percent Error Calc
    P_N{n}.err = abs(P_N{n}.mat - P_true)./(P_true) * 100;
    PerrVec(n) = norm(P_true - P_N{n}.mat, 'fro')./norm(P_true, 'fro')...
                 * 100;
    VerrVec(n) = norm(V_true - V_N{n}.mat, 'fro')./norm(V_true, 'fro')...
                 * 100;
    nVec(n) = N;
end

%% Plotting error contour 
% Some locations incur greater errors than others, and I'd like to show
% how that varies

set(0, 'defaulttextinterpreter', 'latex');
n = 2;
if mod(n-1, 2) == 0 % Deciding number of vortices
    N = 5 * 10^((n-1)/2);
else
    N = 10^((n)/2);
end

% Contour levels
Pmin = min(min(P_N{n}.err)); Pmax = max(max(P_N{n}.err));
levsP = linspace(Pmin, Pmax, 15);

figure
hold on
[C, l] = contourf(P_N{n}.xMat, P_N{n}.yMat, P_N{n}.err, levsP);
set(l,'LineColor','none')
plot([0 c], [0 0], 'r', 'linewidth', 4) % the "airfoil"
title(sprintf('Percent Error for N = %.e Vortices',  nVec(n)))
ylabel('y Position')
xlabel('x Position')
axis([-c/4 c/4 -c/4 c/4]) % Zooming in on the goodies
set(gca, 'TickLabelInterpreter', 'latex')

cBar = colorbar;
cBar.Label.String = 'Percent Error';
cBar.Label.Interpreter = 'latex';
cBar.Label.FontSize = 12;
cBar.TickLabelInterpreter = 'latex';
colormap jet
hold off

%% Error trend

% Fits for error
P_obj = fit(nVec, PerrVec, 'power2');
V_obj = fit(nVec, VerrVec, 'power2');
fitVec = linspace(nVec(1)*0.6, nVec(end)*2, 10000);

figure
datV = semilogx(nVec ,VerrVec, 'o', 'LineWidth', 1, ...
                                    'Color', [0.8, 0, 0]);
hold on % Must go after first semilogx call to format axis correctly
fitV = semilogx(fitVec, V_obj(fitVec), '--r','LineWidth', 1.5);
datP = semilogx(nVec, PerrVec, 'o', 'LineWidth', 1, ...
                                    'Color', [0, 0, 0.8]);
fitP = semilogx(fitVec, P_obj(fitVec), '--b', 'LineWidth', 1.5);


% Figure formatting
leg = legend('Velocity Field Percent Error',...
             'Velocity Error Fit',...
             'Pressure Field Percent Error', ...
             'Pressure Error Fit');
leg.Interpreter = 'latex';
title('Percent Error vs. Number of Vortices')
xlabel('Number of Vortices')
ylabel('Percent Error')
set(gca, 'TickLabelInterpreter', 'latex',...
         'FontSize', 12);
axis([nVec(1)*0.6 nVec(end)*2 0 max(max(PerrVec), max(VerrVec))*1.1]);
grid on
grid minor
hold off
% End figure formatting

toc
fprintf('Velocity Error Trend fit: E = %.1f * N^%.1s + %.1s\n',...
    V_obj.a, V_obj.b, V_obj.c)
fprintf('Pressure Error Trend fit: E = %.1f * N^%.1s + %.1s\n',...
    P_obj.a, P_obj.b, P_obj.c)


