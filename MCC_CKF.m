close all
clear
clc
%% call system dynamic
sysinfo;

%% initiallization
tf = 300;
iter = floor(tf/T); % length of signal
num_of_exprement = 100 ; % number of itterations
num_shot_noise = 15;
start_of_shotnoise = 60;
index_rand_shot = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];

%% Choice the type of noise
SE_KF= zeros(num_of_exprement,num_vec,iter);
SE_CKF_UMV =zeros(num_of_exprement,num_vec,iter);
SE_MCC_CKF=zeros(num_of_exprement,num_vec,iter);
SE_CEnKF =zeros(num_of_exprement,num_vec,iter);
SE_CUKF =zeros(num_of_exprement,num_vec,iter);

for Numexper = 1:num_of_exprement
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
    
    Q = Q_n1;  R = R_n1;
    MeasErrX = sqrt(Q)*randn(num_vec,iter);
    MeasErrZ = sqrt(R)*randn(num_meas,iter);
    
    % initial conditions for system and estimator(x0 , xhat0)
    x = initial_x; % real states
    xhat1 = initial_x; % (standard KF)
    xhat2 = initial_x; % (CKF_UMV)
    xhat3 = initial_x; % (MCC_CKF)
    xhat4 = initial_x; %  (CEnKF)
    xhat5 = initial_x; %  (CUKF)
        
    % elements of initial estimation error covariance (P0hat)
    P_KF = initial_P;  %  (standard KF)
    P_CKF_UMV = initial_P;  %  (CKF_UMV)
    P_MCC_CKF = initial_P; % (MCC_CKF)
    P_CEnKF = initial_P;  %  (CEnKF)
    P_CUKF = initial_P;  %  (CUKF)
    
    %% define some arrays to store the main signal and the esitimation signals
    xhat_KF = zeros(num_vec,iter); % (standard KF)
    xhat_CKF_UMV = zeros(num_vec,iter); %(CKF_UMV)
    xhat_MCC_CKF = zeros(num_vec,iter); % (MCC_CKF)
    xhat_CEnKF = zeros(num_vec,iter); % (CEnKF)
    xhat_CUKF = zeros(num_vec,iter); % (CUKF)
    
    
    x_main = zeros(num_vec,iter);
    %     x_main(:,1) = x ; %  (main signal)
    
    Shot_Noise;
    
     % number of samples in Enkf
    N_Enkf = 100;
    no = sqrt(P_CEnKF) * repmat(randn(size(xhat4)),1,N_Enkf);
    X_CEnKF = repmat(xhat4,1,N_Enkf) + no;
    sigma = 5;
    for t = 1 : 1 : iter
        % make the measurement signal
        z = B*x;
        z = z + MeasErrZ(:,t);
        
        %%  ======================= KF ========================
        xhat1=A*xhat1;
        P_KF=A*P_KF*A'+Q;
        
        K=P_KF*B'*inv(B*P_KF*B'+R);
        xhat1=xhat1+K*(z-B*xhat1);
        P_KF=(eye(4)-K*B)*P_KF;
        
        xhat_KF(:,t) = xhat1;
        
        %%  ======================= CKF-UMV ========================
        K=P_CKF_UMV*B'*inv(B*P_CKF_UMV*B'+R);
        xhat2=xhat2+K*(z-B*xhat2);
        P_CKF_UMV=(eye(4)-K*B)*P_CKF_UMV;
        
        Gain = P_CKF_UMV*D'*inv(D*P_CKF_UMV*D');
        xhat2=xhat2 - Gain*(D*xhat2 - d);
        P_CKF_UMV = (eye(num_vec) - Gain*D)*P_CKF_UMV*(eye(num_vec) - Gain*D)';
        xhat_CKF_UMV(:,t) =  xhat2;
        
        xhat2=A*xhat2;
        P_CKF_UMV=A*P_CKF_UMV*A'+Q;
        
        %% ======================= run MCC_CKF ========================
        xhat3 = A * xhat3;
        P_MCC_CKF = A * P_MCC_CKF  * A' + Q;
        invers_R = pinv(R);
        innov = z - B * xhat3;
        sigma = 1000 ;
        W = 1;
        norm_innov1 = (innov(1))'*invers_R(1,1)*(innov(1));
        C11 = exp(-(norm_innov1^2) /(2*sigma^2));
        norm_innov2 = (innov(2))'*invers_R(2,2)*(innov(2));
        C22 = exp(-(norm_innov2^2) /(2*sigma^2));
        Cm = diag([C11 C22]);
        lambda = pinv(P_MCC_CKF) + B'* Cm* invers_R * B - sigma^2* D' * inv(W) * D;
        L1 = pinv(lambda)*B'* Cm* invers_R;
        L2 = pinv(lambda)*sigma^2* D' * inv(W);
        
        
%%

            
%%
 
%        G1 = diag(exp(-(diag((diag(innov.^2)*invers_R))./(2*sigma^2))));
%        eig(B'*G1*inv(R)*B  - B'*inv(R)*B)
%        sigma = max(diag(abs(sqrt(inv(diag(diag(D'*D))) * (B'*G1*inv(R)*B  - B'*inv(R)*B)))))
%         C_x = diag([1;1;1;1]);
%         P_cx = C_x * P_MCC_CKF;
%         R_cy = Cm * inv(R);
%         Gain = pinv( B' * R_cy * B + pinv(P_MCC_CKF)+sigma^2*D'*D);
%         Gain1 = Gain * B' * R_cy;
%         Gain2 = sigma^2 * Gain * D';
        xhat3 = xhat3 + L1 *(innov) + L2 * (D*xhat3 - d);
        temp = (eye(num_vec) - L1*B + L2 * D);
        P_MCC_CKF  = temp *P_MCC_CKF *temp' + L1 * R *L1';
        xhat_MCC_CKF(:,t)=xhat3;
        
        %% =======================(CEnKF) =====================
        X_CEnKF = A * X_CEnKF;
        no = sqrt(Q) * randn(size(X_CEnKF));
        X_CEnKF = X_CEnKF +  no;
%         for i = 1 : N_Enkf
%    X_CEnKF(:,i) = X_CEnKF(:,i) - P_CEnKF * D' * pinv(D*P_CEnKF*D') * (D * X_CEnKF(:,i) - d);
% end
        xhat4 = mean(X_CEnKF,2);
        P_CEnKF = (1/(N_Enkf-1))*((X_CEnKF - repmat(xhat4,1,N_Enkf))*(X_CEnKF - repmat(xhat4,1,N_Enkf)).');
        K = P_CEnKF*B'/(B*P_CEnKF*B'+R);
        X_CEnKF = X_CEnKF + K * (repmat(z,1,N_Enkf)- B * X_CEnKF);
        xhat4 = mean(X_CEnKF,2);
        P_CEnKF = (1/(N_Enkf-1))*((X_CEnKF - repmat(xhat4,1,N_Enkf))*(X_CEnKF-repmat(xhat4,1,N_Enkf)).');
        xhat4 = xhat4 - P_CEnKF * D' * inv(D*P_CEnKF*D') * (D * xhat4 - d);
        
        xhat_CEnKF(:,t)=xhat4;
        
        %% ======================= (CUKF) ========================
        [xhat5,P_CUKF] = ukf(F,xhat5,P_CUKF,H,z,Q,R,D,d);
        xhat5 = xhat5 - P_CUKF * D' * inv(D*P_CUKF*D') * (D * xhat5 - d);
        xhat_CUKF(:,t)=xhat5;
        %% Simulate the system.
        x = A*x + MeasErrX(:,t);
        
        %% Simulate the system.
        
        x = A * x +  sqrt(Q)*randn(size(x));
        % Constrain the vehicle (i.e., the true state) to the straight road.
        if abs(x(1) - tan(teta) * x(2)) > 5
            x(2) = (x(2) + x(1) * tan(teta)) / (1 + tan(teta)^2);
            x(1) = x(2) * tan(teta);
        end
        if abs(x(3) - tan(teta) * x(4)) > 0.2
            x(4) = (x(4) + x(3) * tan(teta)) / (1 + tan(teta)^2);
            x(3) = x(4) * tan(teta);
        end
               x_main(:,t+1) = x;
        
    end
    x_main(:,iter) = [];
    
    SE_KF(Numexper,:,:)=(x_main - xhat_KF).^2;
    SE_CKF_UMV(Numexper,:,:) =(x_main - xhat_CKF_UMV).^2;
    SE_MCC_CKF(Numexper,:,:)=(x_main - xhat_MCC_CKF).^2;
    SE_CEnKF(Numexper,:,:) =(x_main - xhat_CEnKF).^2;
    SE_CUKF(Numexper,:,:) =(x_main - xhat_CUKF).^2;
    
end
%%
MSE_KF = zeros(num_vec,iter);
MSE_CKF_UMV= zeros(num_vec,iter);
MSE_MCC_CKF= zeros(num_vec,iter);
MSE_CEnKF= zeros(num_vec,iter);
MSE_CUKF= zeros(num_vec,iter);

for i = 1 : iter
    MSE_KF(:,i) = sqrt(mean(SE_KF(:,:,i)))';
    MSE_CKF_UMV(:,i) = sqrt(mean(SE_CKF_UMV(:,:,i)))';
    MSE_MCC_CKF(:,i) = sqrt(mean(SE_MCC_CKF(:,:,i)))';
    MSE_CEnKF(:,i) = sqrt(mean(SE_CEnKF(:,:,i)))';
    MSE_CUKF(:,i) = sqrt(mean(SE_CUKF(:,:,i)))';
    
end
MMSE_KF = mean(MSE_KF,2);
MMSE_MCC_CKF = mean(MSE_MCC_CKF,2);
MMSE_CKF_UMV = mean(MSE_CKF_UMV,2);
MMSE_CEnKF = mean(MSE_CEnKF,2);
MMSE_CUKF = mean(MSE_CUKF,2);
%

MMMSE_KF = mean(MMSE_KF);
MMMSE_CKF_UMV = mean(MMSE_CKF_UMV);
MMMSE_MCC_CKF = mean(MMSE_MCC_CKF);
MMMSE_CEnKF = mean(MMSE_CEnKF);
MMMSE_CUKF = mean(MMSE_CUKF);
%% Normlized RMSE
max_KF= max(MSE_KF.');
max_CKF_UMV= max(MSE_CKF_UMV.');
max_MCC_CKF= max(MSE_MCC_CKF.');
max_CEnKF= max(MSE_CEnKF.');
max_CUKF= max(MSE_CUKF.');
NRMSE = zeros(1,5);

j=1;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_KF(i,1)/ max_KF(1,i);
end
j=2;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_CKF_UMV(i,1)/ max_CKF_UMV(1,i);
end
j = 3;
for i=1:4
    NRMSE(1,j) =NRMSE(1,j) + MMSE_MCC_CKF(i,1)/ max_MCC_CKF(1,i);
end
j=4;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_CEnKF(i,1)/ max_CEnKF(1,i);
end
j=5;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_CUKF(i,1)/ max_CUKF(1,i);
end


%% Plot data.
close all
SetPlotOptions;
t = 0 : T : tf-T;
c = 3;

figure,
xlabel('time'), ylabel('Value of shot noise')

hold on
X = index_rand_shot;
Y = MeasErrZ(1,index_rand_shot(1:end));
for i=1:num_shot_noise
    bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
end

figure
hold on
plot(t,MSE_KF(c,:) ,'black',t,MSE_CKF_UMV(c,:) ,'blue',t,MSE_MCC_CKF(c,:) ,'g');legend('KF','CKF-UMV','MCC-CKF'),xlabel('time'), ylabel('Root mean square error')

%% Table result

disp('Mean square error : ');
disp('           x1          x2          x3        x4');
disp(['KF    : ',num2str(MMSE_KF.'),'']);
disp(['CKF_UMV  : ',num2str(MMSE_CKF_UMV.'),'']);
disp(['MCC_CKF : ',num2str(MMSE_MCC_CKF.'),'']);
disp(['CEnKF   : ',num2str(MMSE_CEnKF.'),'']);
disp(['CUKF   : ',num2str(MMSE_CUKF.'),'']);

