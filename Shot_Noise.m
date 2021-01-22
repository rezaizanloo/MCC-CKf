global MeasErry
%%%% Generate Non-guassian Noise(Shot Noise)
% %% In process equation
%     for i=1:num_shot_noise
%         temp = randi([1 5],num_vec,1);
%         MeasErrX(:,index_rand_shot(i)) =  MeasErrX(:,index_rand_shot(i)) + temp + 0.6*randn(num_vec,1); %#ok<SAGROW>
%     end
%  Q = cov(MeasErrX');

%% In Measurment equation
for i=1:num_shot_noise
    %     temp = randi([1 30],num_meas,1);
    %     MeasErrZ(:,index_rand_shot(i)) =  MeasErrZ(:,index_rand_shot(i)) + temp + 0.6*randn(num_meas,1); %#ok<SAGROW>
    
    temp = randi([1 1000],1,1);
    %% shot noise just for temp2
        temp2 = randi([1 2],1,1);
%     temp2 = 1;
    MeasErry(temp2,index_rand_shot(i)) =  MeasErry(temp2,index_rand_shot(i)) + temp + 0.6*randn(1,1); %#ok<SAGROW>
end
