function GLM_coeff = genRecode(spike_directory,directory,destination)
% function GLM_coeff = genRecode(spike_directory,directory,destination)
%
%
% INPUTS:
%
% OUTPUTS:

cd(spike_directory)

cell_files = dir('*.mat');
cue_mod(1:443) = 0;

for kk = 1:length(dir('*.mat'))
    load(cell_files(kk).name);
    mat_overview.fname{kk} = cell_files(kk).name;
    disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
    
    block_drift.block1_length(kk) = length(FRATE.Cue.Trial_firing_rate_block1);
    block_drift.block1_half(kk) = round(block_drift.block1_length(kk) / 2);
    block_drift.b1_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)));
    block_drift.b1_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    block_drift.MWU_b1(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    
    block_drift.block2_length(kk) = length(FRATE.Cue.Trial_firing_rate_block2);
    block_drift.block2_half(kk) = round(block_drift.block2_length(kk) / 2);
    block_drift.b2_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)));
    block_drift.b2_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    block_drift.MWU_b2(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    
    switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
        case 0
            if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                cue_mod(kk) = 1 ;
            end
    end
end

%%
cd(directory)

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome
Predictors = {'Modality' 'Location' 'Outcome'};
for iEpoch = 1:length(Epoch)
    mat_files = dir('GLM_*');
    iWindow = -.5:.1:.5;
    for iGLM = 1:length(iWindow)
        current_file = strcat('GLM_',Epoch{iEpoch},'_DATA_',num2str(iWindow(iGLM)),'.mat');
        for iFind = 1:length(mat_files)
            if strcmp(mat_files(iFind).name,current_file) == 1
                load(strcat(directory,mat_files(iFind).name),'ALL_matrix');
                for iPred = 1:length(Predictors)
                    disp(cat(2,'Epoch ',num2str(iEpoch),' (file #',num2str(iGLM),') Pred ',num2str(iPred)));
                    count = 1;
                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(1:133,iGLM) = 0;
                    for iCMod = 1:length(cue_mod)
                        if iCMod <= length(ALL_matrix)
                            if cue_mod(iCMod) == 1
                                if ALL_matrix(iCMod,iPred) == 1
                                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(count,iGLM) = 1;
                                end
                                count = count + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%%
for iEpoch = 1:3
    for iEpoch2 = 1:3
        for iGLM = 1:11
            for iGLM2 = 1:11
                for iPred = 1:3
                    [GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues] = corrcoef(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2));
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr(2);
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues(2);
                end
            end
        end
    end
end

%%
for iEpoch = 1:3
    for iGLM = 1:11
        for iGLM2 = 1:11
            for iPred = 1:3
                for iPred2 = 1:3
                    [GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues] = corrcoef(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2})(:,iGLM2));
                    GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr(2);
                    GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues(2);
                end
                
            end
        end
    end
end

%%
start = 6; %use 1 or 6, when you want start of analysis

for iEpoch = 1:3
    
    GLM_coeff.summary.(Epoch{iEpoch}).Corr = [];
    GLM_coeff.summary.(Epoch{iEpoch}).Pvalue = [];
    for iPred = 1:3
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue = [];
        if iEpoch == 1
            GLM_coeff.summary.(Predictors{iPred}).Corr = [];
            GLM_coeff.summary.(Predictors{iPred}).Pvalue = [];
        end
        for iEpoch2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(start:11,start:11));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(start:11,start:11));
        end
        GLM_coeff.summary.(Predictors{iPred}).Corr = cat(1,GLM_coeff.summary.(Predictors{iPred}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr);
        GLM_coeff.summary.(Predictors{iPred}).Pvalue = cat(1,GLM_coeff.summary.(Predictors{iPred}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue);
        for iPred2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr(start:11,start:11));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues(start:11,start:11));
        end
        GLM_coeff.summary.(Epoch{iEpoch}).Corr = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr);
        GLM_coeff.summary.(Epoch{iEpoch}).Pvalue = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue);
    end
end

save(cat(2,destination,'Correlation_matrices_DATA.mat'),'GLM_coeff','cue_mod','GLM_recode')

end