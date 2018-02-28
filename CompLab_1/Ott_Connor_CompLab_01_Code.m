%{
Both Problems for Comp Lab 1 are computed in this script. The print
statement at the end holds all my answers. 
%}


%% Part 1 Code
%{
Performs the calculations for Problem 1 of Comp Lab 1

Connor Ott
Created: 9/1/17
Last Modified: 9/14/17
%}

clear variables
close all
clc

% Determining Cp over the cylinder
Cp = @(x) 1 - 4*(sin(x).^2);
xVec = linspace(0, 2*pi, 100); % Just for looking at the plot
CpVec = Cp(xVec);

% Determining Pressure Distribution
rho_inf = 0.9093; % kg/m^3
P_inf = 7.012e4; % Pa
V_inf = 25; % m/s

P_dist = @(x) (Cp(x)* 0.5 * rho_inf * V_inf^2) + P_inf;

% Getting lift and drag based on cylinder geometry
Lift_func = @(x) -P_dist(x) .* sin(x);
Drag_func = @(x) -P_dist(x) .* cos(x);

% Beginning Numerical Integration
N = 4; % (To be changed at a later date)
h = 2 * pi / N;

cumLift = 0; % Cumulative lift and drag
cumDrag = 0;
% Summation using Simpson's rule. 
for k = 1:N/2
    t2km1 = (2*k - 2) * h;
    t2kp1 = (2*k) * h;
    t2k = (2*k - 1) * h;
    
    cumLift = cumLift + Lift_func(t2km1) + 4*(Lift_func(t2k)) + Lift_func(t2kp1);
    cumDrag = cumDrag + Drag_func(t2km1) + 4*(Drag_func(t2k)) + Drag_func(t2kp1);
end

% Lift and Drag Calc.
lift = cumLift * h / 3  * 2;
drag = cumDrag * h / 3 * 2;

fprintf('Lift on the Cylinder: %.3s N\n', lift);
fprintf('Drag on the Cylinder: %.3s N\n\n', drag);
fprintf('Calculated using %.f intervals\n\n', N);

%% Part 2 Code
%{
Performs calculations for part 2 of Comp Lab 1

Created: 9/4/17 - Connor Ott
Last Modified: 9/14/17 - Connor Ott
%}

clear variables
close all

alpha = 9 * pi/180; % angle of attack [rad]

%loading in splines and changing them from Cp to P.
load cp
Cp_lower.coefs = Cp_lower.coefs * 0.5 * 1.225 * 20^2 + 10.13e4;
Cp_upper.coefs = Cp_upper.coefs * 0.5 * 1.225 * 20^2 + 10.13e4;

%% Creating the airfoil
t = 12/100; % maximum thickness as fraction of chord
c = 0.5; % Chord length.
shapeFunc = @(x) t/0.2 * c * (0.2969*(x/c).^(1/2) - 0.1260 * (x/c) - ...
    0.3516*(x/c).^2 + 0.2843*(x/c).^3 - 0.1036*(x/c).^4);
xVecFunc = linspace(0, c, 1000);
shapeVec = shapeFunc(xVecFunc);

%% Summation - Finding Actual Lift and Drag
N_act = 5e5; % # subintervals to Obtain actual Values. (500,000)
xVec = linspace(0, c, N_act+1);
shapeVec = shapeFunc(xVec);
p_vecUp = fnval(Cp_upper, xVec/c); % Pressure vectors for upper and lower surfaces.
p_vecLo = fnval(Cp_lower, xVec/c);


normalSum_lo = 0;
axialSum_lo = 0;
normalSum_hi = 0;
axialSum_hi = 0;

for k = 1:N_act
    deltaX = xVec(k+1) - xVec(k);
    xk = xVec(k); % x at k
    xkp1 = xVec(k + 1); % x at k plus 1
    deltaY = shapeVec(k + 1) - shapeVec(k);
    
    % Normal Forces
    normalSum_hi = normalSum_hi - (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaX;
    normalSum_lo = normalSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaX;
    
    % Axial Forces
    axialSum_hi = axialSum_hi + (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaY;
    axialSum_lo = axialSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaY;
end

% Summing upper and lower surface
normalSum = normalSum_hi + normalSum_lo;
axialSum = axialSum_hi + axialSum_lo;

% Determining Lift and Drag
lift_act = normalSum * cos(alpha) - axialSum * sin(alpha);
drag_act = normalSum * sin(alpha) + axialSum * cos(alpha);

%% Determining subintervals required for different accuracy tolerances.

for i = 1:3
    N = 20; % Starting N for each case.
    lift = 0; % Starting Lift for while loop to evaluate.
    switch i % Switching between tolerances
        case 1 % case for tol < 5%.
            while (lift_act - lift)/lift_act >= 0.05
                xVec = linspace(0, c, N+1);
                shapeVec = shapeFunc(xVec);
                p_vecUp = fnval(Cp_upper, xVec/c);
                p_vecLo = fnval(Cp_lower, xVec/c);
                
                normalSum_lo = 0;
                axialSum_lo = 0;
                normalSum_hi = 0;
                axialSum_hi = 0;
                
                for k = 1:N
                    deltaX = xVec(k+1) - xVec(k);
                    xk = xVec(k); % x at k
                    xkp1 = xVec(k + 1); % x at k plus 1
                    deltaY = shapeVec(k + 1) - shapeVec(k);
                    
                    % Normal Forces
                    normalSum_hi = normalSum_hi - (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaX;
                    normalSum_lo = normalSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaX;
                    
                    % Axial Forces
                    axialSum_hi = axialSum_hi + (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaY;
                    axialSum_lo = axialSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaY;
                end
                % Summing upper and lower surface
                normalSum = normalSum_hi + normalSum_lo;
                axialSum = axialSum_hi + axialSum_lo;
                
                % Determining Lift and Drag
                lift = normalSum * cos(alpha) - axialSum * sin(alpha);
                drag = normalSum * sin(alpha) + axialSum * cos(alpha);
                N = N + 1;
            end
            N_5per = N - 1; % minus one to account for that last plus 1
        case 2 % Case for tol < 1.
            while (lift_act - lift)/lift_act >= 0.01
                xVec = linspace(0, c, N+1);
                shapeVec = shapeFunc(xVec);
                p_vecUp = fnval(Cp_upper, xVec/c);
                p_vecLo = fnval(Cp_lower, xVec/c);
                
                normalSum_lo = 0;
                axialSum_lo = 0;
                normalSum_hi = 0;
                axialSum_hi = 0;
                
                for k = 1:N
                    deltaX = xVec(k+1) - xVec(k);
                    xk = xVec(k); % x at k
                    xkp1 = xVec(k + 1); % x at k plus 1
                    deltaY = shapeVec(k + 1) - shapeVec(k);
                    
                    % Normal Forces
                    normalSum_hi = normalSum_hi - (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaX;
                    normalSum_lo = normalSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaX;
                    
                    % Axial Forces
                    axialSum_hi = axialSum_hi + (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaY;
                    axialSum_lo = axialSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaY;
                end
                % Summing upper and lower surface
                normalSum = normalSum_hi + normalSum_lo;
                axialSum = axialSum_hi + axialSum_lo;
                
                % Determining Lift and Drag
                lift = normalSum * cos(alpha) - axialSum * sin(alpha);
                drag = normalSum * sin(alpha) + axialSum * cos(alpha);
                N = N + 1; % trying the next N
            end
            N_1per = N - 1; % minus one to account for that last plus 1
        case 3 % case for tol < 0.1%
            while (lift_act - lift)/lift_act >= 0.001
                xVec = linspace(0, c, N+1);
                shapeVec = shapeFunc(xVec);
                p_vecUp = fnval(Cp_upper, xVec/c);
                p_vecLo = fnval(Cp_lower, xVec/c);
                
                normalSum_lo = 0;
                axialSum_lo = 0;
                normalSum_hi = 0;
                axialSum_hi = 0;
                
                for k = 1:N
                    deltaX = xVec(k+1) - xVec(k);
                    xk = xVec(k); % x at k
                    xkp1 = xVec(k + 1); % x at k plus 1
                    deltaY = shapeVec(k + 1) - shapeVec(k);
                    
                    % Normal Forces
                    normalSum_hi = normalSum_hi - (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaX;
                    normalSum_lo = normalSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaX;
                    
                    % Axial Forces
                    axialSum_hi = axialSum_hi + (p_vecUp(k) + p_vecUp(k+1))*0.5*deltaY;
                    axialSum_lo = axialSum_lo + (p_vecLo(k) + p_vecLo(k+1))*0.5*deltaY;
                end
                % Summing upper and lower surface
                normalSum = normalSum_hi + normalSum_lo;
                axialSum = axialSum_hi + axialSum_lo;
                
                % Determining Lift and Drag
                lift = normalSum * cos(alpha) - axialSum * sin(alpha);
                drag = normalSum * sin(alpha) + axialSum * cos(alpha);
                N = N + 1;
            end
            N_p1per = N - 1; % point 1 percent - % minus one to account for that last plus 1
        otherwise
            % do nothing you useless piece of trash.
    end
    
end

fprintf('\nLift calculated from %.s intervals on NACA 0012: %.2f N\n', N_act, lift_act);
fprintf('Drag calculated from %.s intervals on NACA 0012: %.2f N\n\n', N_act, drag_act);
fprintf('Subintervals required for 5%% accuracy: %.f\n', N_5per);
fprintf('Subintervals required for 1%% accuracy: %.f\n', N_1per);
fprintf('Subintervals required for 0.1%% accuracy: %.f\n', N_p1per);




