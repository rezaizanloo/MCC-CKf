function [EstErrors, ConstrErrors] = Pendulum(nMonte, DisplayFlag, CorrectQFlag)

% This function is used to study constrained Kalman filtering for the
% nonlinear pendulum system with nonlinear constraints.
% INPUTS:   nMonte = # of Monte Carlo simulations to run
%           DisplayFlag = flag indicating whether or not to display results during run
%           CorrectQFlag = flag indicating whether or not to use the corrected Q (via system projection)
% OUTPUTS:  See function NonlinConstr for explanation of outputs

global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0 MeasErry

if ~exist('nMonte', 'var')
    nMonte = 1;
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('CorrectQFlag', 'var')
    CorrectQFlag = false;
end

tf = 5; % final time (seconds)
T = 0.05; % time step (seconds)
L = 1; % pendulum length
g = 9.81; % acceleration due to gravity
m = 1; % pendulum mass
n = 2; % number of states
num_meas = 2;
System = @PendulumSystem;
Constraint = @PendulumConstraint;

horizon = 2; % MHE horizon length - set this to 0 to skip the MHE code
%% 
iter = floor(tf/T); % length of signal
num_shot_noise = 15;
start_of_shotnoise = 2;
index_rand_shot = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];

% Measurement matrix
H = [1 0; 0 1];

Q = diag([0.007^2, 0.007^2]); % Process noise covariance (m, m, m/sec, m/sec)
R = diag([0.1^2, 0.1^2]); % Measurement noise covariance (m, m)

% Define the initial state and the initial estimation error covariance
x = [pi/4; pi/50];
P0 = eye(2);
y = H * x;
y = y + sqrt(R) * randn(size(y));
E0 = Constraint(x);

EstErrors = zeros(nMonte, 11);
ConstrErrors = zeros(nMonte, 11);
for i = 1 : nMonte
    disp(['Simulation # ', num2str(i)]);
   MeasErry = sqrt(R)*randn(num_meas,iter);
   Shot_Noise;
    [EstErrors(i,:), ConstrErrors(i,:)] = ...
        NonlinConstr(tf, x, P0, @Jacobian, @ConstraintDeriv, DisplayFlag, CorrectQFlag);
end
disp('RMS Estimation Errors(unconstr, perfect meas, est proj, sys proj, pdf trunc, SCKF, UKF, Proj UKF, Eq UKF, 2-step UKF,MCC-KF):');
disp(mean(EstErrors, 1));
% disp('RMS Constraint Errors (unconstr, perfect meas, est proj, MHE, sys proj, pdf trunc, SCKF, UKF, Proj UKF, Eq UKF, 2-step UKF,MCC-KF):');
% disp(mean(ConstrErrors, 1));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, E] = PendulumSystem(x)
% Integration of one time step for pendulum system
global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0
% Trapezoidal integration
dx1 = [x(2); -g / L * sin(x(1))] * T;
dx2 = [x(2) + dx1(2); -g / L * sin(x(1)+dx1(1))] * T;
x = x + (dx1 + dx2) / 2;
y = H * x;
E = Constraint(x);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E] = PendulumConstraint(x)
% Pendulum constraint equation
global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0
E = m * L * L / 2 * x(2)^2 - m * g * L * cos(x(1));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F] = Jacobian(x)
% Jacobian for pendulum system equation
global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0
F = [1 T;
    -T*g/L*cos(x(1)) 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E1, E2, d] = ConstraintDeriv(x)
% Derivatives of constraint for pendulum system
global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0
% E1 is the first derivative of the constraint w/r to the states.
% This is the D matrix obtained via linearization.
E1 = [m*g*L*sin(x(1)), m*L*L*x(2)]; 
% E2 is the second derivative of the constraint w/r to the states.
E2 = [m*g*L*cos(x(1)), 0; 
    0, m*L*L];
% d is the constrained value of Dx = d, obtained via linearization.
d = E0 - Constraint(x) + E1 * x;
return;
