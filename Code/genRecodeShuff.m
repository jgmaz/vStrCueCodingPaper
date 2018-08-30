function GLM_coeff_SHUFF = genRecode(directory,destination,num_Shuffs)
% function GLM_coeff_SHUFF = genRecode(directory,destination,num_Shuffs)
%
%
% INPUTS:
%
% OUTPUTS:

load(cat(2,directory,'Correlation_matrices_DATA.mat'))
cd(directory)

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome
Predictors = {'Modality' 'Location' 'Outcome'};

%% across epochs
rng('shuffle')
num_GLM = 11; %6 or 11 depending if want all time windows

for iEpoch = 1:3
    for iEpoch2 = 1:3
        for iGLM = 1:num_GLM
            for iGLM2 = 1:num_GLM
                disp(cat(2,'Epochs ',num2str(iEpoch),' (',num2str(iGLM),') & ',num2str(iEpoch2),' (',num2str(iGLM2),')'))
                for iPred = 1:3
                    for iShuff = 1:num_Shuffs
                        Epoch1_shuff = [];
                        Epoch2_shuff = [];
                        Epoch1_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM)),'Replace',false);
                        Epoch2_shuff = datasample(GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2),length(GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2)),'Replace',false);
                        [GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(Epoch1_shuff,Epoch2_shuff);
                        GLM_coeff_SHUFF.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr(2);
                        GLM_coeff_SHUFF.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(2);
                    end
                end
            end
        end
    end
end

%% across predictors
for iEpoch = 1:3
    for iGLM = 1:num_GLM
        for iGLM2 = 1:num_GLM
            for iPred = 1:3
                disp(cat(2,'Epoch ',num2str(iEpoch),' (',num2str(iGLM),') & (',num2str(iGLM2),')'))
                for iPred2 = 1:3
                    for iShuff = 1:num_Shuffs
                        Pred1_shuff = [];
                        Pred2_shuff = [];
                        Pred1_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM)),'Replace',false);
                        Pred2_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2})(:,iGLM2),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2})(:,iGLM2)),'Replace',false);
                        [GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(Pred1_shuff,Pred2_shuff);
                        GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr(2);
                        GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(2);
                    end
                end
            end
        end
    end
end

%%
start = 6; %use 1 or 6, when you want start of analysis

for iShuff = 1:num_Shuffs
    disp(iShuff)
    for iEpoch = 1:3
        
        GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr = [];
        GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue = [];
        for iPred = 1:3
            GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr = [];
            GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue = [];
            GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr = [];
            GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue = [];
            if iEpoch == 1
                GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr = [];
                GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue = [];
            end
            for iEpoch2 = 1:3
                GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr = cat(2,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorr(start:num_GLM,start:num_GLM));
                GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue = cat(2,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(start:num_GLM,start:num_GLM));
            end
            GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr = cat(1,GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr);
            GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue = cat(1,GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue);
            for iPred2 = 1:3
                GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr = cat(2,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorr(start:num_GLM,start:num_GLM));
                GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue = cat(2,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(start:num_GLM,start:num_GLM));
            end
            GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr = cat(1,GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr);
            GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue = cat(1,GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue,GLM_coeff_SHUFF.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue);
        end
    end
end

%% Z-scores away from shuffle (value - shuff mean / shuff std)
%%%% for MEAN of Corr %%%%
matrix_start = [1 7 13];
for iShuff = 1:num_Shuffs
    for iEpoch = 1:3
        for iRow = 1:3
            for iCol = 1:3
                tempEpoch = GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr(matrix_start(iRow):matrix_start(iRow)+5,matrix_start(iCol):matrix_start(iCol)+5);
                Table.Shuffs.(Epoch{iEpoch}){iRow,iCol}(iShuff) = mean(tempEpoch(:));
            end
        end
    end
    
    for iPred = 1:3
        for iRow = 1:3
            for iCol = 1:3
                tempPred = GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr(matrix_start(iRow):matrix_start(iRow)+5,matrix_start(iCol):matrix_start(iCol)+5);
                Table.Shuffs.(Predictors{iPred}){iRow,iCol}(iShuff) = mean(tempPred(:));
            end
        end
    end
end

for iEpoch = 1:3
    for iRow = 1:3
        for iCol = 1:3
            tempEpoch = GLM_coeff.summary.(Epoch{iEpoch}).Corr(matrix_start(iRow):matrix_start(iRow)+5,matrix_start(iCol):matrix_start(iCol)+5);
            Table.(Epoch{iEpoch}).Data(iRow,iCol) = mean(tempEpoch(:));
            Table.(Epoch{iEpoch}).ShuffMEAN(iRow,iCol) = mean(Table.Shuffs.(Epoch{iEpoch}){iRow,iCol});
            Table.(Epoch{iEpoch}).ShuffSTD(iRow,iCol) = std(Table.Shuffs.(Epoch{iEpoch}){iRow,iCol});
            Table.(Epoch{iEpoch}).Zscore(iRow,iCol) = (Table.(Epoch{iEpoch}).Data(iRow,iCol) - Table.(Epoch{iEpoch}).ShuffMEAN(iRow,iCol)) / Table.(Epoch{iEpoch}).ShuffSTD(iRow,iCol);
            if Table.(Epoch{iEpoch}).Zscore(iRow,iCol) > 1.96
                Table.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = 1;
            elseif Table.(Epoch{iEpoch}).Zscore(iRow,iCol) < -1.96
                Table.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = -1;
            else
                Table.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = 0;
            end
        end
    end
end

for iPred = 1:3
    for iRow = 1:3
        for iCol = 1:3
            tempPred = GLM_coeff.summary.(Predictors{iPred}).Corr(matrix_start(iRow):matrix_start(iRow)+5,matrix_start(iCol):matrix_start(iCol)+5);
            Table.(Predictors{iPred}).Data(iRow,iCol) = mean(tempPred(:));
            Table.(Predictors{iPred}).ShuffMEAN(iRow,iCol) = mean(Table.Shuffs.(Predictors{iPred}){iRow,iCol});
            Table.(Predictors{iPred}).ShuffSTD(iRow,iCol) = std(Table.Shuffs.(Predictors{iPred}){iRow,iCol});
            Table.(Predictors{iPred}).Zscore(iRow,iCol) = (Table.(Predictors{iPred}).Data(iRow,iCol) - Table.(Predictors{iPred}).ShuffMEAN(iRow,iCol)) / Table.(Predictors{iPred}).ShuffSTD(iRow,iCol);
            if Table.(Predictors{iPred}).Zscore(iRow,iCol) > 1.96
                Table.(Predictors{iPred}).Zscore_recode(iRow,iCol) = 1;
            elseif Table.(Predictors{iPred}).Zscore(iRow,iCol) < -1.96
                Table.(Predictors{iPred}).Zscore_recode(iRow,iCol) = -1;
            else
                Table.(Predictors{iPred}).Zscore_recode(iRow,iCol) = 0;
            end
        end
    end
end

save(cat(2,destination,'Correlation_matrices_SHUFF.mat'),'GLM_coeff_SHUFF','Table')

end