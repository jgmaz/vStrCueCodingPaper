%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
count_Value = 1;
for iCell = 1:133
GLM_allTime{iCell}(1:54) = 0;
end
% Predictors = {'Modality' 'Location' 'Outcome'};
for iEpoch = 1:3%:length(Epoch)
    switch iEpoch
        case 1
            mat_files = dir('2018-03-24-GLM_cueon*');
            Predictors = {'Modality' 'Location' 'Outcome'}; % 'Approach' 'Latency' 'Trial' 'Previous'};
        case 2
            mat_files = dir('2018-03-26-GLM_NP*');
            Predictors = {'Modality' 'Location' 'Outcome'};
        case 3
            mat_files = dir('2018-03-26-GLM_outcome*');
            Predictors = {'Modality' 'Location' 'Outcome'};
%         case 4
%             mat_files = dir('2018-03-12-GLM_cueoff*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-03-12-GLM_cueoff_-0.5-500bin.mat' '2018-03-12-GLM_cueoff_-0.4-500bin.mat' '2018-03-12-GLM_cueoff_-0.3-500bin.mat' ...
%     '2018-03-12-GLM_cueoff_-0.2-500bin.mat' '2018-03-12-GLM_cueoff_-0.1-500bin.mat' '2018-03-12-GLM_cueoff_0-500bin.mat' '2018-03-12-GLM_cueoff_0.1-500bin.mat' ...
%     '2018-03-12-GLM_cueoff_0.2-500bin.mat' '2018-03-12-GLM_cueoff_0.3-500bin.mat' '2018-03-12-GLM_cueoff_0.4-500bin.mat' '2018-03-12-GLM_cueoff_0.5-500bin.mat'};
    end
    % for iGLM = 1:length(file_order)
    iWindow = -.5:.1:.5;
    for iGLM = 6:length(iWindow)
        if iEpoch == 1
            current_file = cat(2,'2018-03-24-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'.mat');
%         elseif iEpoch == 4
%             current_file = cat(2,'2018-03-12-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'-500bin.mat');
        else
            current_file = cat(2,'2018-03-26-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'-round1.mat');
        end
        for iFind = 1:length(mat_files)
            if strcmp(mat_files(iFind).name,current_file) == 1
                load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files(iFind).name),'ALL_matrix');
                
                for iPred = 1:length(Predictors)
                    disp(cat(2,'Epoch ',num2str(iEpoch),' (file #',num2str(iGLM),') Pred ',num2str(iPred)));
                    count_Cell = 1;
                    for iCMod = 1:length(cue_mod)
                        if iCMod <= length(ALL_matrix)
                            if cue_mod(iCMod) == 1
                                if ALL_matrix(iCMod,iPred) == 1
                                    GLM_allTime{count_Cell}(count_Value) = 1;                                    
                                end
                                count_Cell = count_Cell + 1;
                            end
                        end
                    end
                    count_Value = count_Value + 1;
                end
            end
        end
    end
end

%%
cue_presence(1:133) = 0;
for iCell = 1:133
    GLM_test(iCell,1) = sum(GLM_allTime{iCell});
    if GLM_test(iCell,1) > 0
        cue_presence(iCell) = 1;
    end
end