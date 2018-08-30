function GLM_window = genGLMstats(destination,num_Shuffs,select_Epoch)
% function GLM_window = genGLMstats(destination,num_Shuffs,select_Epoch)
%
%
% INPUTS:
%
% OUTPUTS:

cd(destination)
start_point = [1 2 1];
end_point = [1 3 3];

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome
DataType = {'DATA' 'SHUFF'};
for iEpoch = start_point(select_Epoch):end_point(select_Epoch)
    disp(iEpoch)
    find_files.DATA = strcat('GLM_',Epoch{iEpoch},'_DATA');
    find_files.SHUFF = strcat('GLM_',Epoch{iEpoch},'_SHUFF');
    if iEpoch == 1
        predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
        predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
    end
    
    mat_files.DATA = dir(strcat(find_files.DATA,'*'));
    mat_files.SHUFF = dir(strcat(find_files.SHUFF,'*'));
    iTime = -0.5:.1:0.5;
    for iFile = 1:length(mat_files.DATA)
        for iType = 1:2
            file_order.(DataType{iType}){iFile} = strcat(find_files.(DataType{iType}),'_',num2str(iTime(iFile)),'.mat');
        end
    end
    
    for iGLM = 1:length(file_order.DATA)
        for iFind = 1:length(mat_files.DATA)
            if strcmp(mat_files.DATA(iFind).name,file_order.DATA{iGLM}) == 1
                load(mat_files.DATA(iFind).name);
                disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files.DATA))));
                temp_Rsquared = GLM_matrices.Rsquared.ALL;
                temp_Rsquared(temp_Rsquared == 0) = NaN;
                temp_Rsquared(temp_Rsquared > 50) = NaN;
                temp_Rsquared(temp_Rsquared < -50) = NaN;
                temp_Rsquared = abs(temp_Rsquared);
                for iList = 1:length(predictor_list)
                    GLM_window.DATA.(Epoch{iEpoch}).window(iList,iGLM) = summary_var.All.(predictor_list{iList});
                    if iList > 1 && iList < 12
                        GLM_window.DATA.(Epoch{iEpoch}).Rsquared.MEAN(iList-1,iGLM) = nanmean(temp_Rsquared(:,iList-1));
                        GLM_window.DATA.(Epoch{iEpoch}).Rsquared.SEM(iList-1,iGLM) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                    end
                end
                
            end
        end
    end
    
    GLM_window.DATA.(Epoch{iEpoch}).prop = GLM_window.DATA.(Epoch{iEpoch}).window / 133;
    
    for iGLM = 1:length(file_order.SHUFF)
        for iFind = 1:length(mat_files.SHUFF)
            if strcmp(mat_files.SHUFF(iFind).name,file_order.SHUFF{iGLM}) == 1
                load(mat_files.SHUFF(iFind).name,'GLM_matrices','summary_var');
                disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files.SHUFF))));
                for iShuff = 1:num_Shuffs
                    temp_Rsquared = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL;
                    temp_Rsquared(temp_Rsquared == 0) = NaN;
                    temp_Rsquared(temp_Rsquared > 50) = NaN;
                    temp_Rsquared(temp_Rsquared < -50) = NaN;
                    temp_Rsquared = abs(temp_Rsquared);
                    for iList = 1:length(predictor_list)
                        GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).window(iList,iGLM) = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.(predictor_list{iList});
                        if iList > 1 && iList < 12
                            GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.MEAN(iList-1,iGLM) = nanmean(temp_Rsquared(:,iList-1));
                            GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.SEM(iList-1,iGLM) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                        end
                    end
                end
            end
        end
    end
    
    % end
    
    %%
    % for iEpoch = start_point(select_Epoch):end_point(select_Epoch)
    for iShuff = 1:num_Shuffs
        GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).prop = GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).window / 133;
        for iList = 2:length(predictor_list)
            GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})(iShuff,:) = GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).prop(iList,:);
            GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList})(iShuff,:) = GLM_window.SHUFF.ALL.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.MEAN(iList-1,:);
        end
    end
    for iList = 2:length(predictor_list)
        GLM_window.SHUFF.(Epoch{iEpoch}).prop.MEAN.(predictor_list{iList}) = nanmean(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList}));
        GLM_window.SHUFF.(Epoch{iEpoch}).prop.SEM.(predictor_list{iList}) = nanstd(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList}))/ ...
            sqrt(numel(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})(:,1)));
        GLM_window.SHUFF.(Epoch{iEpoch}).prop.ZSTD.(predictor_list{iList}) = nanstd(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})) *1.96;
        GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.MEAN.(predictor_list{iList}) = nanmean(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList}));
        GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.SEM.(predictor_list{iList}) = nanstd(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList}))/ ...
            sqrt(numel(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList})(:,1)));
    end
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Z-scores away from shuffle (value - shuff mean / shuff std)
    % for iEpoch = start_point(select_Epoch):end_point(select_Epoch);
    if iEpoch == 1
        Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
        Predictors = {'Modality' 'Location' 'Outcome'};
    end
    
    for iPred = 1:length(Predictors);
        GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Data = GLM_window.DATA.(Epoch{iEpoch}).prop(iPred+1,:);
        GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).ShuffMEAN = mean(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(Predictors{iPred}));
        GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).ShuffSTD = std(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(Predictors{iPred}));
        for iTime = 1:length(GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Data)
            GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Zscore(iTime) = (GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Data(iTime) - GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).ShuffMEAN(iTime)) / GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).ShuffSTD(iTime);
            if GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Zscore(iTime) > 1.96
                GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Zscore_recode(iTime) = 1;
            else
                GLM_window.Table.(Epoch{iEpoch}).(Predictors{iPred}).Zscore_recode(iTime) = 0;
            end
        end
    end
end

switch select_Epoch
    case 1
        save(cat(2,destination,'GLM_cueon_plotting_stats.mat'),'GLM_window')
    case 2
        save(cat(2,destination,'GLM_NP-out_plotting_stats.mat'),'GLM_window')
    case 3
        save(cat(2,destination,'GLM_allEpochs_plotting_stats.mat'),'GLM_window')
end

end