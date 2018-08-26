% cell_files = dir('*.mat');
count = 1;
for iCMod = 1:length(cue_mod)
    
    if cue_mod(iCMod) == 1
        
        for iShuff = 1:100
            if iCMod <= length(ALL_matrix.(cat(2,'shuff_',num2str(iShuff))))
                WithinCorr{count}(iShuff,1) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,1);
                WithinCorr{count}(iShuff,2) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,2);
                WithinCorr{count}(iShuff,3) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,3);
                
                
                
            end
            
        end
        count = count + 1;
    end
end

%%

% for iCell = 1:length(WithinCorr)
%     WithCorr{iCell} = corrcoef(WithinCorr{iCell});
%     Corrs(iCell,1) = WithCorr{iCell}(2);
%     Corrs(iCell,2) = WithCorr{iCell}(3);
%     Corrs(iCell,3) = WithCorr{iCell}(6);
%     
% end

num_Shuffs = 100;
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
% Predictors = {'Modality' 'Location' 'Outcome'};
for iEpoch = 1:3%:length(Epoch)
    switch iEpoch
        case 1
            mat_files = dir('2018-03-17-GLM_cueon*');
            Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
        case 2
            mat_files = dir('2018-03-17-GLM_NP*');
            Predictors = {'Modality' 'Location' 'Outcome'};
        case 3
            mat_files = dir('2018-03-17-GLM_outcome*');
            Predictors = {'Modality' 'Location' 'Outcome'};
    end
    % for iGLM = 1:length(file_order)
    iWindow = -.5:.1:.5;
    for iGLM = 1:length(iWindow)
        current_file = cat(2,'2018-03-17-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'-SHUFF.mat');
        for iFind = 1:length(mat_files)
            if strcmp(mat_files(iFind).name,current_file) == 1
                load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files(iFind).name),'ALL_matrix');
                disp(cat(2,'Epoch ',num2str(iEpoch),' (file #',num2str(iGLM),')'));
                count = 1;
                for iCMod = 1:length(cue_mod)
                    if cue_mod(iCMod) == 1
                        for iShuff = 1:100
                            if iCMod <= length(ALL_matrix.(cat(2,'shuff_',num2str(iShuff))))
                                WithinCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){count}(iShuff,1) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,1);
                                WithinCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){count}(iShuff,2) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,2);
                                WithinCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){count}(iShuff,3) = ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,3);
                            end
                        end
                        count = count + 1;
                    end
                end
            end
        end
    end
end

%%
for iEpoch = 1:3
    for iGLM = 1:length(iWindow)

for iCell = 1:length(WithinCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))))
    WithCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){iCell} = corrcoef(WithinCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){iCell});
    Corrs.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM)))(iCell,1) = WithCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){iCell}(2);
    Corrs.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM)))(iCell,2) = WithCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){iCell}(3);
    Corrs.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM)))(iCell,3) = WithCorr.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))){iCell}(6);
    
end
MeanCorr.(Epoch{iEpoch})(iGLM,:) = nanmean(Corrs.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))));
STDCorr.(Epoch{iEpoch})(iGLM,:) = nanstd(Corrs.(Epoch{iEpoch}).(cat(2,'w',num2str(iGLM))));
    end
end