%% examples  information

%% Example 1 info
teta = pi/3;
T = 3 ;
num_vec = 4; % dimension of state vector
num_meas = 2; % number of measurment
A =[1 0 T 0  % process Equation
    0 1 0 T
    0 0 1 0
    0 0 0 1];
B =[1 0 0 0  % meaurment equation
    0 1 0 0];
D = [1 -tan(teta) 0 0;0 0 1 -tan(teta)];
d = [0;0];
initial_x = [1;1;0;0]; % initial of state vector X
initial_P = diag([4,4,3,3]);% initial of cov matrix P
mu_n1_x = [-3;-3;-3;-3]; % mean of noise1 state vector
mu_n2_x = [2;2;2;2]; % mean of noise2 state vector
mu_n1_z = [-2;-2]; % mean of noise1 measurment
mu_n2_z = [2;2]; % mean of noise2 measurment
R_n1 = diag([0.1,0.1]);% cov of Noise1 measurment
R_n2 = diag([0.1,0.1]); % cov of Noise2 measurment
Q_n1 = diag([0.0001,0.0001,0.0001,0.0001]); % cov of Noise1 process
Q_n2 = diag([0.0001,0.0001,0.0001,0.0001]); % cov of Noise2 process
F = @(x)[x(1)+T*x(3);x(2)+T*x(4);x(3);x(4)];  % nonlinear state equations for UKF
H = @(x)[x(1);x(2)];                               % measurement equation for UKF

%% Example 2 info
% teta = pi/3;
% global T
% T = 3 ;
% num_vec = 6; % dimension of state vector
% num_meas = 2; % number of measurment
% A=[1 T T^2/2 0 0 0 % process equation
%     0 1 T 0 0 0
%     0 0 1 0 0 0
%     0 0 0 1 T T^2/2
%     0 0 0 0 1 T
%     0 0 0 0 0 1];
% B =[1 0 0 0 0 0  % measurment equation
%     0 0 0 1 0 0];
% initial_x = [1;1;0;0;1;1]; % initial of state vector X
% initial_P = diag([4,4,3,3,4,4]);% initial of cov matrix P
% mu_n1_x = [-2;-2;-2;-2;-2;-2]; % mean of noise1 state vector
% mu_n2_x = [0;0;0;0;0;0]; % mean of noise2 state vector
% mu_n1_z = [0;0]; % mean of noise1 measurment
% mu_n2_z = [1;1]; % mean of noise2 measurment
% R_n1 = diag([0.01,0.01]);% cov of Noise1 measurment
% R_n2 = diag([0.01,0.01]); % cov of Noise2 measurment
% Q_n1 = diag([0.01,0.01,0.01,0.01,0.01,0.01]); % cov of Noise1 process
% Q_n2 = diag([0.01,0.01,0.01,0.01,0.01,0.01]); % cov of Noise2 process
% F = @(x)[x(1)+x(2)+x(3)/2;x(2)+T*x(3);x(3);x(4)+x(5)+x(6)/2;x(5)+x(6);x(6)];  % nonlinear state equations for UKF
% H = @(x)[x(1);x(4)];                               % measurement equation for UKF

%% Example 3 info
% teta = pi/3;
% num_vec = 2; % dimension of state vector
% num_meas = 1; % number of measurment
% A=[cos(teta) sin(teta);   % process Equation   [0.5000    0.8660
%     -sin(teta) cos(teta)];    %                 -0.8660    0.5000]
% B =[1 1];  % measurment equation
% initial_x = [1;1]; % initial of state vector X
% initial_P = diag([4,4]);% initial of cov matrix P
% mu_n1_x = [0;0]; % mean of noise1 state vector
% mu_n2_x = [0;0]; % mean of noise2 state vector
% mu_n1_z = [0]; % mean of noise1 measurment
% mu_n2_z = [0]; % mean of noise2 measurment
% Q_n1 = [0.01 0 ;0 0.01]; % cov of Noise1 measurment
% Q_n2 = [0.01 0 ;0 0.01 ]; % cov of Noise2 measurment
% R_n1 = [0.01 ]; % cov of Noise1 process
% R_n2 = [0.01 ]; % cov of Noise2 process
% F = @(x)[cos(teta)*x(1)+sin(teta)*x(2);-sin(teta)*x(1)+cos(teta)*x(2)];  % nonlinear state equations for UKF
% H = @(x)[x(1)+x(2)];                               % measurement equation for UKF
