num_Shuffs = 100;
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
 Predictors = {'Modality' 'Location' 'Outcome'};

count = 1;
for iEpoch = 1:3
    for iPred = 1:length(Predictors)
        disp(iPred)
        for iShuff = 1:num_Shuffs
            for iCell = 1:133
                GLM_sum{iCell}(count) = sum(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(iCell,:));
                GLM_sum_exp{iCell}(count,:) = GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(iCell,:);
            end
            count = count + 1;
        end
    end
end

%%
% GLM_include(1:133) = 0;
% for iCell = 1:133
%     GLM_counts(iCell,1) = sum(GLM_sum{iCell});
%     if GLM_counts(iCell) > 990
%         GLM_include(iCell) = 1;
%     end
% end

%%
% for iCell = 1:133
%     GLM_extended(iCell,1) = sum(GLM_sum{iCell}(1:300));
%     GLM_extended(iCell,2) = sum(GLM_sum{iCell}(301:600));
%     GLM_extended(iCell,3) = sum(GLM_sum{iCell}(601:900));
%     GLM_extended(iCell,4) = sum(GLM_sum{iCell}([1:100 301:400 601:700]));
%     GLM_extended(iCell,5) = sum(GLM_sum{iCell}([101:200 401:500 701:900]));
%     GLM_extended(iCell,6) = sum(GLM_sum{iCell}([201:300 501:600 801:900]));
% end

%%
for iCell = 1:133
    temp = sum(GLM_sum_exp{iCell}(1:300,6:11));
    temp = cat(2,temp,sum(GLM_sum_exp{iCell}(301:600,6:11)));
    GLM_exp_extended(iCell,:) = cat(2,temp,sum(GLM_sum_exp{iCell}(601:900,6:11)));
end

corrs_epoch = corrcoef(GLM_exp_extended);

%%
matrix_start = [1 7 13];
    for iRow = 1:3
        for iCol = 1:3
            tempEpoch = corrs_epoch(matrix_start(iRow):matrix_start(iRow)+5,matrix_start(iCol):matrix_start(iCol)+5);
            corrs_Epoch(iRow,iCol) = mean(tempEpoch(:));
        end
    end
%%
% for iCell = 1:133
%     temp = sum(GLM_sum_exp{iCell}([1:100 301:400 601:700],:));
%     temp = cat(2,temp,sum(GLM_sum_exp{iCell}([1:100 301:400 601:700],:)));
%     GLM_exp_extended_cue(iCell,:) = cat(2,temp,sum(GLM_sum_exp{iCell}([1:100 301:400 601:700],:)));
% end
% 
% corrs_pred = corrcoef(GLM_exp_extended_cue);

GLM_means(:,1) = mean(GLM_exp_extended(:,1:6),2);
GLM_means(:,2) = mean(GLM_exp_extended(:,7:12),2);
GLM_means(:,3) = mean(GLM_exp_extended(:,13:18),2);

%%
GLM_exclude(1:133) = 1;
for iCell = 1:133
    temp = GLM_means(iCell,1) / GLM_means(iCell,2);
%     if temp > 3 || temp < .33 || temp == 0 || temp == Inf
        if temp > .67 && temp < 1.5
        GLM_exclude(iCell) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% across epochs within units
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'};
 Predictors = {'Modality' 'Location' 'Outcome'};

 
 GLM_recode_shuff.RAW = [];
 
for iEpoch = 1:3
    for iPred = 1:3
        GLM_recode_shuff.RAW = cat(2,GLM_recode_shuff.RAW,GLM_recode_DATA.(Epoch{iEpoch}).(Predictors{iPred})(:,6:11));
    end
end

%%
rng('shuffle')
num_Shuffs = 100;

for iCell = 1:133
    disp(iCell)
    for iShuff = 1:num_Shuffs
        GLM_recode_shuff.(cat(2,'s',num2str(iShuff)))(iCell,:) = datasample(GLM_recode_shuff.RAW(iCell,:),length(GLM_recode_shuff.RAW(iCell,:)),'Replace',false);
    end
end

%%
col_start = 1;
for iEpoch = 1:3
    for iPred = 1:3
        disp(iPred)
        for iShuff = 1:num_Shuffs
            GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))) = GLM_recode_shuff.(cat(2,'s',num2str(iShuff)))(:,col_start:col_start+5);
        end
        col_start = col_start + 6;
    end
end