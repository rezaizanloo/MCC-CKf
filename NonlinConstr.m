function [EstErrors, ConstrErrors] = ...
    NonlinConstr(tf, x, P0, Jacobian, ConstraintDeriv, DisplayFlag, CorrectQFlag)

% function NonlinConstr
% This m-file simulates nonlinear Kalman filtering with nonlinear state constraints.
% INPUTS:   tf = simulation time
%           horizon = MHE horizon length
%           x = initial state
%           P0 = covariance of initial estimation error
%           Jacobian = address of function to compute Jacobian of system
%           ConstraintDeriv = address of function to compute derivative of state constraint
%           DisplayFlag = flag indicating whether or not to display results during run
%           CorrectQFlag = flag indicating whether or not to use corrected Q (via system projection)
% OUTPUTS:  EstErrors = Array of RMS estimation errors. This is an array containing results for:
%               (1) Unconstrained Kalman filter
%               (2) Perfect measurement filter
%               (3) Estimate projection filter (W=P^{-1})
%               (4) System projection filter
%               (5) Pdf truncation filter
%               (6) Smoothly constrained filter
%               (7) Unscented Kalman filter
%               (8) Projected UKF
%               (9) Equality constrained UKF
%               (10) Two-step UKF
%               (11) MCC-KF
%           ConstrErrors = Array of RMS constraint errors.

global m g T L E0 n H Q R System Constraint xhatMHE0 yMHEArray PMHE0 MeasErry

Rsqrt = sqrt(R);
xhat = x; % Unconstrained Kalman filter
xhat1 = x; % Kalman filter with perfect measurements
xtilde = x; % Constrained Kalman filter (estimate projection, W=inv(P))
xhatSys = x; % Constrained Kalman filter (system projection)
xTrunc = x; % Constrained Kalman filter (pdf truncation)
xhatSCKF = x; % Smoothly constrained Kalman filter
xhatUKF = x; % Unscented Kalman filter
xhatUKFProj = x; % Projected UKF
xhatUKFEq = x; % Equality constrained UKF
xhatUKF2 = x; % Two-step constrained UKF
xhatMCC = x; % MCC_CKF

% Measurement noise covariance for perfect measurement filter
R1 = R;
for i = 1 : length(E0)
    R1(end+1, end+1) = 0;
end

% Initial estimation error covariance for constrained Kalman filter (system projection)
[PND] = NullProj(xhatSys, ConstraintDeriv);
Pc = PND * P0 * PND;
% Process noise covariance for constrained Kalman filter (system projection).
Qc = PND * Q * PND;

if CorrectQFlag
    Q = Qc; % Use the correct P0 and Q for all the filters
    P0 = Pc;
end

[dQc, lambdaQc, dQcT] = svd(Qc);
for i = 1 : size(lambdaQc, 1)
    if real(lambdaQc(i,i)) < 0
        % This should never be true, but it may occur because of numerical issues.
        lambdaQc(i,i) = 0;
    end
end
ddT = dQc * dQc';
if (norm(eye(size(ddT)) - ddT) > 1e-8)
    disp('Error - dQc is not orthogonal.');
    return;
end
if (norm(dQc*lambdaQc*dQc' - Qc) > 1e-8)
    disp('Error - SVD failed for Qc');
    return;
end

% Initial estimation error covariances
P = P0;
P1 = P0;
PUKF = P0;
PUKFEq = P0;
P_MCC_CKF = P0;

% SCKF parameters
sThreshold = 100;
alpha = 0.01;
PSCKF = P0;

% Initialize arrays for saving data for plotting.
xarray = x;
xhatarray = xhat;
xhat1array = xhat1;
xtildearray = xtilde;
xhatSysarray = xhatSys;
xTruncArray = xTrunc;
xhatSCKFArray = xhatSCKF;
xhatUKFArray = xhatUKF;
xhatUKFProjArray = xhatUKFProj;
xhatUKFEqArray = xhatUKFEq;
xhatUKF2Array = xhatUKF2;
xhatMCCArray = xhatMCC;
ConstrXArray = Constraint(x);
Constrarray = Constraint(xhat) - E0;
Constr1array = Constraint(xhat1) - E0;
ConstrTildearray = Constraint(xtilde) - E0;
ConstrSysarray = Constraint(xhatSys) - E0;
ConstrTruncArray = Constraint(xTrunc) - E0;
ConstrSCKFArray = Constraint(xhatSCKF) - E0;
ConstrUKFArray = Constraint(xhatUKF) - E0;
ConstrUKFProjArray = Constraint(xhatUKFProj) - E0;
ConstrUKFEqArray = Constraint(xhatUKFEq) - E0;
ConstrUKF2Array = Constraint(xhatUKF2) - E0;
ConstrMCCArray = Constraint(xhatMCC) - E0;

randn('state', sum(100*clock));
In = eye(n, n);

% Begin the simulation.
for t = T : T : tf

    %% Simulate the system.
    [x, y, E] = System(x);
    y = y + MeasErry(:,floor(t/T));
    %% Run the Kalman filter.
    F = Jacobian(xhat);
    P = F * P * F' + Q;
    xhat = System(xhat);
    K = P * H' * inv(H * P * H' + R);
    xhat = xhat + K * (y - H * xhat);
    P = (In - K * H) * P;
    %% Find the constrained (estimate projection) Kalman filter estimates.
    [D, E2, d] = ConstraintDeriv(xhat);
    xtilde = xhat - P * D' * inv(D*P*D') * (D * xhat - d);
    %% Run the constrained Kalman filter (system projection).
    F = Jacobian(xhatSys);
    PND = NullProj(xhatSys, ConstraintDeriv);
    Qc = PND * Q * PND;
    Pc = F * Pc * F' + Qc;
    [xhatSys] = System(xhatSys);
    Kc = Pc * H' * inv(H * Pc * H' + R);
    xhatSys = xhatSys + Kc * (y - H * xhatSys);
    Pc = (In - Kc * H) * Pc;
    %% Run the filter for the perfect measurement formulation.
    F = Jacobian(xhat1);
    [D, E2, d] = ConstraintDeriv(xhat1);
    y1 = [y; d];
    P1 = F * P1 * F' + Q;
    xhat1 = System(xhat1);
    H1 = [H; D];
    K1 = P1 * H1' * inv(H1 * P1 * H1' + R1);
    xhat1 = xhat1 + K1 * (y1 - H1 * xhat1);
    P1 = (In - K1 * H1) * P1;
   
    %% Constrained filtering via probability density function truncation
    PTrunc = P;
    xTrunc = xhat;
    [D, E2, d] = ConstraintDeriv(xTrunc);    
    for k = 1 : length(E0)
        [Utrunc, Wtrunc, Vtrunc] = svd(PTrunc);
        Ttrunc = Utrunc;
        TTT = Ttrunc * Ttrunc';
        if (norm(eye(size(TTT)) - TTT) > 1e-8)
           disp('Error - Ttrunc is not orthogonal.');
           return;
        end
        if (norm(Utrunc*Wtrunc*Utrunc' - PTrunc) > 1e-8)
            disp('Error - SVD failed for pdf trunction');
            return;
        end
        % Compute the modified Gram-Schmidt transformation S * Amgs = [ Wmgs ; 0 ].
        % Amgs is a given n x m matrix, and S is an orthogonal n x n matrix, and Wmgs is an m x m matrix.
        Amgs = sqrt(Wtrunc) * Ttrunc' * D(k,:)'; % n x 1, where n = number of states
        [Wmgs, S] = MGS(Amgs);
        S = S * sqrt(D(k,:) * PTrunc * D(k,:)') / Wmgs;
        cTrunc = (d(k) - D(k,:) * xTrunc) / sqrt(D(k,:) * PTrunc * D(k,:)');
        dTrunc = (d(k) - D(k,:) * xTrunc) / sqrt(D(k,:) * PTrunc * D(k,:)');
        % The next 3 lines are commented out because they apply only for inequality constraints
        %alpha = sqrt(2/pi) / erf(dTrunc/sqrt(2)) - erf(cTrunc/sqrt(2));
        %mu = alpha * (exp(-cTrunc^2/2) - exp(-dTrunc^2/2));
        %sigma2 = alpha * (exp(-cTrunc^2/2) * (cTrunc - 2 * mu) - exp(-dTrunc^2/2) * (dTrunc - 2 * mu)) + mu^2 + 1;
        
        % Since we have equality constraints, the E(zTrunc) = cTrunc = dTrunc,
        % and var(zTrunc) = 0.
        mu = cTrunc;
        sigma2 = 0;
        
        zTrunc = zeros(size(xTrunc));
        zTrunc(1) = mu;
        CovZ = eye(length(zTrunc));
        CovZ(1,1) = sigma2;
        xTrunc = Ttrunc * sqrt(Wtrunc) * S' * zTrunc + xTrunc;
        PTrunc = Ttrunc * sqrt(Wtrunc) * S' * CovZ * S * sqrt(Wtrunc) * Ttrunc';
    end
    %% Smoothly constrained Kalman filter
    F = Jacobian(xhatSCKF);
    PSCKF = F * PSCKF * F' + Q;
    xhatSCKF = System(xhatSCKF);
    KSCKF = PSCKF * H' * inv(H * PSCKF * H' + R);
    xhatSCKF = xhatSCKF + KSCKF * (y - H * xhatSCKF);
    PSCKF = (In - KSCKF * H) * PSCKF;
    E1 = ConstraintDeriv(xhatSCKF);
    zeta0 = alpha * E1 * PSCKF * E1';
    for NumApp = 0 : 10
        zeta = zeta0 * exp(-NumApp);
        sSCKF = max(E1' .* diag(PSCKF) .* E1') / (E1 * PSCKF * E1');
        if (sSCKF > sThreshold), break, end
        KSCKF = PSCKF * E1' * inv(E1 * PSCKF * E1' + zeta);
        xhatSCKF = xhatSCKF + KSCKF * (E0 - Constraint(xhatSCKF));
        PSCKF = (In - KSCKF * E1) * PSCKF * (In - KSCKF * E1)' + KSCKF * zeta * KSCKF';        
        E1 = ConstraintDeriv(xhatSCKF);
    end
    %% Unscented Kalman filters
    [xhatUKF, PUKF, xhatUKFProj, temp1, xhatUKF2] = UKF(xhatUKF, PUKF, System, Constraint, ConstraintDeriv, E0, H, Q, R, y);
    [temp1, temp2, xhatUKFEq, PUKFEq] = UKF(xhatUKFEq, PUKFEq, System, Constraint, ConstraintDeriv, E0, H, Q, R, y);
    %% MCC_CKF
        xhatMCC = System(xhatMCC);
        F = Jacobian(xhatMCC);
        P_MCC_CKF = F * P_MCC_CKF  * F' + Q;
         E1 = ConstraintDeriv(xhatMCC);
        invers_R = pinv(R);
        innov = y - H * xhatMCC;
        sigma = 5 ;
       G1 = diag(exp(-(diag((diag(innov.^2)*invers_R))./(2*sigma^2))));
        R_cy = G1 * inv(R);
        R_cy_2 = G1 * inv(R) * G1';
        Gain = pinv(H' * R_cy * H + pinv(P_MCC_CKF)+sigma^2*E1'*E1);
        Gain1 = Gain * H' * (R_cy);
        Gain2 = sigma^2 * Gain * E1';
        xhatMCC = xhatMCC + Gain1 *(innov) - Gain2 * (E1*xhatMCC - E0);
        temp = (eye(n) - Gain1*H - Gain2 * E1);
        P_MCC_CKF  = temp *P_MCC_CKF *temp' + Gain1 * R *Gain1';
    
    %% Save data in arrays
    xarray = [xarray x];
    xhatarray = [xhatarray xhat];
    xtildearray = [xtildearray xtilde];
    xhat1array = [xhat1array xhat1];
    xhatSysarray = [xhatSysarray xhatSys]; 
    xTruncArray = [xTruncArray xTrunc];
    xhatSCKFArray = [xhatSCKFArray xhatSCKF];
    xhatUKFArray = [xhatUKFArray xhatUKF];
    xhatUKFProjArray = [xhatUKFProjArray xhatUKFProj];
    xhatUKFEqArray = [xhatUKFEqArray xhatUKFEq];
    xhatUKF2Array = [xhatUKF2Array xhatUKF2];
    xhatMCCArray = [xhatMCCArray xhatMCC];
    ConstrXArray = [ConstrXArray Constraint(x)];
    Constrarray = [Constrarray Constraint(xhat)-E0]; 
    ConstrTildearray = [ConstrTildearray Constraint(xtilde)-E0];
    Constr1array = [Constr1array Constraint(xhat1)-E0];
    ConstrSysarray = [ConstrSysarray Constraint(xhatSys)-E0];
    ConstrTruncArray = [ConstrTruncArray Constraint(xTrunc)-E0];
    ConstrSCKFArray = [ConstrSCKFArray Constraint(xhatSCKF)-E0];
    ConstrUKFArray = [ConstrUKFArray Constraint(xhatUKF)-E0];
    ConstrUKFProjArray = [ConstrUKFProjArray Constraint(xhatUKFProj)-E0];
    ConstrUKFEqArray = [ConstrUKFEqArray Constraint(xhatUKFEq)-E0];
    ConstrUKF2Array = [ConstrUKF2Array Constraint(xhatUKF2)-E0];
      ConstrMCCArray = [ConstrMCCArray Constraint(xhatMCC)-E0];
end
% Compute RMS estimation errors
EstError = xarray - xhatarray;
EstError = mean(EstError(1,:).^2 + EstError(2,:).^2);
EstError = sqrt(EstError);
EstError1 = xarray - xhat1array;
EstError1 = mean(EstError1(1,:).^2 + EstError1(2,:).^2);
EstError1 = sqrt(EstError1);
EstErrorTilde = xarray - xtildearray;
EstErrorTilde = mean(EstErrorTilde(1,:).^2 + EstErrorTilde(2,:).^2);
EstErrorTilde = sqrt(EstErrorTilde);
EstErrorSys = xarray - xhatSysarray;
EstErrorSys = mean(EstErrorSys(1,:).^2 + EstErrorSys(2,:).^2);
EstErrorSys = sqrt(EstErrorSys);
EstErrorTrunc = xarray - xTruncArray;
EstErrorTrunc = mean(EstErrorTrunc(1,:).^2 + EstErrorTrunc(2,:).^2);
EstErrorTrunc = sqrt(EstErrorTrunc);
EstErrorSCKF = xarray - xhatSCKFArray;
EstErrorSCKF = mean(EstErrorSCKF(1,:).^2 + EstErrorSCKF(2,:).^2);
EstErrorSCKF = sqrt(EstErrorSCKF);
EstErrorUKF = xarray - xhatUKFArray;
EstErrorUKF = mean(EstErrorUKF(1,:).^2 + EstErrorUKF(2,:).^2);
EstErrorUKF = sqrt(EstErrorUKF);
EstErrorUKFProj = xarray - xhatUKFProjArray;
EstErrorUKFProj = mean(EstErrorUKFProj(1,:).^2 + EstErrorUKFProj(2,:).^2);
EstErrorUKFProj = sqrt(EstErrorUKFProj);
EstErrorUKFEq = xarray - xhatUKFEqArray;
EstErrorUKFEq = mean(EstErrorUKFEq(1,:).^2 + EstErrorUKFEq(2,:).^2);
EstErrorUKFEq = sqrt(EstErrorUKFEq);
EstErrorUKF2 = xarray - xhatUKF2Array;
EstErrorUKF2 = mean(EstErrorUKF2(1,:).^2 + EstErrorUKF2(2,:).^2);
EstErrorUKF2 = sqrt(EstErrorUKF2);
EstErrorMCC = xarray - xhatMCCArray;
EstErrorMCC = mean(EstErrorMCC(1,:).^2 + EstErrorMCC(2,:).^2);
EstErrorMCC = sqrt(EstErrorMCC);
% Compute constraint errors
r = length(d); % number of constraints
Constr = 0; for i = 1 : r, Constr = Constr + Constrarray(i,:).^2; end
Constr = sqrt(mean(Constr));
Constr1 = 0; for i = 1 : r, Constr1 = Constr1 + Constr1array(i,:).^2; end
Constr1 = sqrt(mean(Constr1));
ConstrTilde = 0; for i = 1 : r, ConstrTilde = ConstrTilde + ConstrTildearray(i,:).^2; end
ConstrTilde = sqrt(mean(ConstrTilde));
ConstrSys = 0; for i = 1 : r, ConstrSys = ConstrSys + ConstrSysarray(i,:).^2; end
ConstrSys = sqrt(mean(ConstrSys));
ConstrTrunc = 0; for i = 1 : r, ConstrTrunc = ConstrTrunc + ConstrTruncArray(i,:).^2; end
ConstrTrunc = sqrt(mean(ConstrTrunc));
ConstrSCKF = 0; for i = 1 : r, ConstrSCKF = ConstrSCKF + ConstrSCKFArray(i,:).^2; end
ConstrSCKF = sqrt(mean(ConstrSCKF));
ConstrUKF = 0; for i = 1 : r, ConstrUKF = ConstrUKF + ConstrUKFArray(i,:).^2; end
ConstrUKF = sqrt(mean(ConstrUKF));
ConstrUKFProj = 0; for i = 1 : r, ConstrUKFProj = ConstrUKFProj + ConstrUKFProjArray(i,:).^2; end
ConstrUKFProj = sqrt(mean(ConstrUKFProj));
ConstrUKFEq = 0; for i = 1 : r, ConstrUKFEq = ConstrUKFEq + ConstrUKFEqArray(i,:).^2; end
ConstrUKFEq = sqrt(mean(ConstrUKFEq));
ConstrUKF2 = 0; for i = 1 : r, ConstrUKF2 = ConstrUKF2 + ConstrUKF2Array(i,:).^2; end
ConstrUKF2 = sqrt(mean(ConstrUKF2));
ConstrMCC = 0; for i = 1 : r, ConstrMCC = ConstrMCC + ConstrMCCArray(i,:).^2; end
ConstrMCC = sqrt(mean(ConstrMCC));
EstErrors = [EstError, EstError1, EstErrorTilde, EstErrorSys, EstErrorTrunc , EstErrorSCKF, EstErrorUKF, EstErrorUKFProj, EstErrorUKFEq, EstErrorUKF2,EstErrorMCC];
ConstrErrors = [Constr, Constr1, ConstrTilde, ConstrSys, ConstrTrunc, ConstrSCKF, ConstrUKF, ConstrUKFProj, ConstrUKFEq, ConstrUKF2,ConstrMCC];
if DisplayFlag
    close all; t = 0 : T : tf;
    figure;
    plot(t, xarray(1, :), ':', t, xarray(2, :), '-');
    title('True State');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, ConstrXArray);
    xlabel('seconds'); ylabel('energy');
    figure;
    plot(t, xarray(1, :) - xhatarray(1, :), ':', t, xarray(2, :) - xhatarray(2, :), '-');
    title('Estimation Error (Unconstrained)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhat1array(1, :), ':', t, xarray(2, :) - xhat1array(2, :), '-');
    title('Estimation Error (Perfect Measurements)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xtildearray(1, :), ':', t, xarray(2, :) - xtildearray(2, :), '-');
    title('Estimation Error (Estimate Projection, W=inv(P))');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatSysarray(1, :), ':', t, xarray(2, :) - xhatSysarray(2, :), '-');
    title('Estimation Error (System Projection)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xTruncArray(1, :), ':', t, xarray(2, :) - xTruncArray(2, :), '-');
    title('Estimation Error (pdf Truncation)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatSCKFArray(1, :), ':', t, xarray(2, :) - xhatSCKFArray(2, :), '-');
    title('Estimation Error (SCKF)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatUKFArray(1, :), ':', t, xarray(2, :) - xhatUKFArray(2, :), '-');
    title('Estimation Error (UKF)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatUKFProjArray(1, :), ':', t, xarray(2, :) - xhatUKFProjArray(2, :), '-');
    title('Estimation Error (Projected UKF)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatUKFEqArray(1, :), ':', t, xarray(2, :) - xhatUKFEqArray(2, :), '-');
    title('Estimation Error (Equality constrained UKF)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
    figure;
    plot(t, xarray(1, :) - xhatUKF2Array(1, :), ':', t, xarray(2, :) - xhatUKF2Array(2, :), '-');
    title('Estimation Error (Two-step UKF)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
   figure;
    plot(t, xarray(1, :) - xhatMCCArray(1, :), ':', t, xarray(2, :) - xhatMCCArray(2, :), '-');
    title('Estimation Error (MCC-Kf)');
    xlabel('seconds'); ylabel('meters'); legend('pos', 'vel');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PND] = NullProj(x, ConstraintDeriv)

[D] = ConstraintDeriv(x);
[u, s, v] = svd(D');
r = size(D,1); % number of constraints
u2 = u(:, r+1:end);
PND = u2 * u2';
return;