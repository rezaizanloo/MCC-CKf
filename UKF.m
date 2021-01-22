function [xhatplus, Pplus, xhatUKFProj, PUKFProj, xhat2] = ...
    UKF(xhatplus, Pplus, SystemSim, Constraint, ConstraintDeriv, d, H, Q, R, y)

% Unscented Kalman filter.
% INPUTS:   xhatplus = a posteriori state estimate at the previous time step
%           Pplus = a posteriori state estimate covariance at the previous time step
%           SystemSim = address of system simulation routine
%           Constraint = address of routine that calculates constraint equation
%           ConstraintDeriv = address of routine that calculates constraint derivatives
%           d = constraint value - g(x)=d
%           H = measurement matrix (linear measurement equation)
%           Q = process noise covariance
%           R = measurement noise covariance
%           y = measurement vector
% OUTPUTS:  xhatplus = a posteriori state estimate at the present time step
%           Pplus = a posteriori state estimate covariance
%           xhatUKFProj = projected UKF estimate
%           PUKFProj = state estimate covariance of projected UKF estimate
%           xhat2 = two-step constrained UKF estimate

% UKF time update
N = length(xhatplus);
r = length(y);
[root,p] = chol(N*Pplus);
for i = 1 : N
    sigma(:,i) = xhatplus + root(i,:)';
    sigma(:,i+N) = xhatplus - root(i,:)';
end
for i = 1 : 2*N
    x = SystemSim(sigma(:,i));
    sigma(:,i) = x;
end
xhatminus = zeros(N,1);
for i = 1 : 2*N
    xhatminus = xhatminus + sigma(:,i) / 2 / N;
end
Pminus = zeros(N,N);
for i = 1 : 2*N
    Pminus = Pminus + (sigma(:,i) - xhatminus) * (sigma(:,i) - xhatminus)' / 2 / N;
end
Pminus = Pminus + Q;
% UKF measurement update
[root,p] = chol(N*Pminus);
for i = 1 : N
    sigma(:,i) = xhatminus + root(i,:)';
    sigma(:,i+N) = xhatminus - root(i,:)';
end
for i = 1 : 2*N
    yukf(:,i) = H * sigma(:,i);
end
yhat = 0;
for i = 1 : 2*N
    yhat = yhat + yukf(:,i) / 2 / N;
end
Py = zeros(r,r);
Pxy = zeros(N,r);
for i = 1 : 2*N
    Py = Py + (yukf(:,i) - yhat) * (yukf(:,i) - yhat)' / 2 / N;
    Pxy = Pxy + (sigma(:,i) - xhatminus) * (yukf(:,i) - yhat)' / 2 / N;
end
Py = Py + R;
Kukf = Pxy * inv(Py);
xhatplus = xhatminus + Kukf * (y - yhat);
Pplus = Pminus - Kukf * Py * Kukf';
% Treat the constraint as a perfect measurement
[root,p] = chol(N*Pplus);
for i = 1 : N
    sigma(:,i) = xhatplus + root(i,:)';
    sigma(:,i+N) = xhatplus - root(i,:)';
end
for i = 1 : 2*N
    D(:,i) = Constraint(sigma(:,i));
end
dhat = 0;
for i = 1 : 2*N
    dhat = dhat + D(:,i) / 2 / N;
end
r = size(D,1); % number of constraints
Pdd = zeros(r,r);
Pxd = zeros(N,r);
for i = 1 : 2*N
    Pdd = Pdd + (D(:,i) - dhat) * (D(:,i) - dhat)' / 2 / N;
    Pxd = Pxd + (sigma(:,i) - xhatplus) * (D(:,i) - dhat)' / 2 / N;
end
Kukf = Pxd * inv(Pdd);
xhatUKFProj = xhatplus + Kukf * (d - dhat);
PUKFProj = Pplus - Kukf * Pdd * Kukf';
% Two-step projection for equality constraints
% First project the pdf onto the constraint
[root,p] = chol(N*Pplus);
for i = 1 : N
    sigma(:,i) = xhatplus + root(i,:)';
    sigma(:,i+N) = xhatplus - root(i,:)';
end
for i = 1 : 2*N
    [D, E2, d] = ConstraintDeriv(sigma(:,i));
    sigma(:,i) = sigma(:,i) - Pplus * D' * inv(D*Pplus*D') * (D * sigma(:,i) - d);
end
% Calculate the mean of the projected pdf
xhat2 = zeros(N,1);
for i = 1 : 2*N
    xhat2 = xhat2 + sigma(:,i) / 2 / N;
end
% Calculate the covariance of the projected pdf
Pdag = zeros(N,N);
for i = 1 : 2*N
    Pdag = Pdag + (sigma(:,i) - xhat2) * (sigma(:,i) - xhat2)' / 2 / N;
end
% Project the mean of the projected pdf onto the constraint
[D, E2, d] = ConstraintDeriv(xhat2);
xhat2 = xhat2 - Pdag * D' * inv(D*Pdag*D') * (D * xhat2 - d);
return
