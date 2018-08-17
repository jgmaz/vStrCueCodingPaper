cell_files = dir('*.mat');
cue_mod(1:443) = 0;
cue_MSN(1:443) = 0;
cue_FSI(1:443) = 0;
cue_FSI_Inc(1:443) = 0;
cue_FSI_Dec(1:443) = 0;
cue_MSN_Inc(1:443) = 0;
cue_MSN_Dec(1:443) = 0;

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
            if RANK.two.Trial > 975 || RANK.two.Trial < 26
                if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                    cue_mod(kk) = 1 ;
                    
                    isi = diff(spk_t);
                    sorted_isi = sort(isi,'descend');
                    if sorted_isi(5) < 2
                        cue_FSI(kk) = 1;
                        %%%%%%% separating inc and dec %%%%%%%
                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                            case 1
                                cue_FSI_Inc(kk) = 1;
                            case 0
                                cue_FSI_Dec(kk) = 1;
                        end
                        %%%%%%%% end %%%%%%%%
                    else
                        cue_MSN(kk) = 1;
                        %%%%%%% separating inc and dec %%%%%%%
                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                            case 1
                                cue_MSN_Inc(kk) = 1;
                            case 0
                                cue_MSN_Dec(kk) = 1;
                        end
                        %%%%%%%% end %%%%%%%%
                    end
                end
            end
    end
end

%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
% Predictors = {'Modality' 'Location' 'Outcome'};
for iEpoch = 1:3%:length(Epoch)
    switch iEpoch
        case 1
            mat_files = dir('2018-03-24-GLM_cueon*');
            Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
        case 2
            mat_files = dir('2018-03-26-GLM_NP*');
            Predictors = {'Modality' 'Location' 'Outcome'};
        case 3
            mat_files = dir('2018-03-26-GLM_outcome*');
            Predictors = {'Modality' 'Location' 'Outcome'};
        case 4
            mat_files = dir('2018-03-12-GLM_cueoff*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-03-12-GLM_cueoff_-0.5-500bin.mat' '2018-03-12-GLM_cueoff_-0.4-500bin.mat' '2018-03-12-GLM_cueoff_-0.3-500bin.mat' ...
%     '2018-03-12-GLM_cueoff_-0.2-500bin.mat' '2018-03-12-GLM_cueoff_-0.1-500bin.mat' '2018-03-12-GLM_cueoff_0-500bin.mat' '2018-03-12-GLM_cueoff_0.1-500bin.mat' ...
%     '2018-03-12-GLM_cueoff_0.2-500bin.mat' '2018-03-12-GLM_cueoff_0.3-500bin.mat' '2018-03-12-GLM_cueoff_0.4-500bin.mat' '2018-03-12-GLM_cueoff_0.5-500bin.mat'};
    end
    % for iGLM = 1:length(file_order)
    iWindow = -.5:.1:.5;
    for iGLM = 1:length(iWindow)
        if iEpoch == 1
            current_file = cat(2,'2018-03-24-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'.mat');
        elseif iEpoch == 4
            current_file = cat(2,'2018-03-12-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'-500bin.mat');
        else
            current_file = cat(2,'2018-03-26-GLM_',Epoch{iEpoch},'_',num2str(iWindow(iGLM)),'-round1.mat');
        end
        for iFind = 1:length(mat_files)
            if strcmp(mat_files(iFind).name,current_file) == 1
                load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files(iFind).name),'ALL_matrix');
                
                for iPred = 1:length(Predictors)
                    disp(cat(2,'Epoch ',num2str(iEpoch),' (file #',num2str(iGLM),') Pred ',num2str(iPred)));
                    count = 1;
%                     count_MSN = 1;
%                     count_FSI = 1;
                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(1:133,iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).ALL(1:sum(cue_mod),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN(1:sum(cue_MSN),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI(1:sum(cue_FSI),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN_Inc(1:sum(cue_MSN_Inc),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI_Inc(1:sum(cue_FSI_Inc),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN_Dec(1:sum(cue_MSN_Dec),iGLM) = 0;
                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI_Dec(1:sum(cue_FSI_Dec),iGLM) = 0;
                    for iCMod = 1:length(cue_mod)
                        if iCMod <= length(ALL_matrix)
                            if cue_mod(iCMod) == 1
                                if ALL_matrix(iCMod,iPred) == 1
                                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(count,iGLM) = 1;
                                    GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).ALL(count,iGLM) = 1;                                    
                                    if cue_MSN(iCMod) == 1
                                        GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN(sum(cue_MSN(1:iCMod)),iGLM) = 1;
                                        if cue_MSN_Dec(iCMod) == 1
                                            GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN_Dec(sum(cue_MSN_Dec(1:iCMod)),iGLM) = 1;
                                        elseif cue_MSN_Inc(iCMod) == 1
                                            GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).MSN_Inc(sum(cue_MSN_Inc(1:iCMod)),iGLM) = 1;
                                        end
%                                         count_MSN = count_MSN + 1;
                                    elseif cue_FSI(iCMod) == 1
                                        GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI(sum(cue_FSI(1:iCMod)),iGLM) = 1;
                                        if cue_FSI_Dec(iCMod) == 1
                                            GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI_Dec(sum(cue_FSI_Dec(1:iCMod)),iGLM) = 1;
                                        elseif cue_FSI_Inc(iCMod) == 1
                                            GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).FSI_Inc(sum(cue_FSI_Inc(1:iCMod)),iGLM) = 1;
                                        end
%                                         count_FSI = count_FSI + 1;
                                    end
                                    
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
Class = {'ALL' 'MSN_Inc' 'MSN_Dec' 'FSI_Inc' 'FSI_Dec'};
cue_type = {'cue_mod' 'cue_MSN_Inc' 'cue_MSN_Dec' 'cue_FSI_Inc' 'cue_FSI_Dec'};
for iEpoch = 1:3
    if iEpoch == 1
         Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
         Predictors = {'Modality' 'Location' 'Outcome'};
    end
            Summary_table.(Epoch{iEpoch})(1:length(Predictors),1:length(Class)) = 0;
    for iPred = 1:length(Predictors)
        for iClass = 1:length(Class)
            Summary_table.(Epoch{iEpoch})(iPred,iClass) = sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,6));
            Summary_percent.(Epoch{iEpoch})(iPred,iClass) = (sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,6))/sum(cue_type{iClass}))*100;
        end
    end
end
%%
Class = {'ALL' 'MSN' 'FSI' 'MSN_Inc' 'MSN_Dec' 'FSI_Inc' 'FSI_Dec'};
Total_count = [sum(cue_mod) sum(cue_MSN) sum(cue_FSI) sum(cue_MSN_Inc) sum(cue_MSN_Dec) sum(cue_FSI_Inc) sum(cue_FSI_Dec)];

for iEpoch = 1:3
    for iEpoch2 = 1:3
        for iGLM = 1:11
            for iGLM2 = 1:11
                for iPred = 1:3
                    [GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues] = corrcoef(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2));
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr(2);
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues(2);
%                     for iClass = 1:length(Class)
%                     GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).CountOverlap(iGLM,iGLM2) = sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM) == 1 & GLM_unitType.RAW.(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass})(:,iGLM2) == 1); 
%                     GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).Proportion(iGLM,iGLM2) = GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).CountOverlap(iGLM,iGLM2) / sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM));
%                     GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).PropTotal(iGLM,iGLM2) = GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).CountOverlap(iGLM,iGLM2) / Total_count(iClass);
%                     
%                     [GLM_unitType.corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr GLM_unitType.corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues] = corrcoef(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM),GLM_unitType.RAW.(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass})(:,iGLM2));
%                     GLM_unitType.coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).RecodeCorr(iGLM,iGLM2) = GLM_unitType.corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr(2);
%                     GLM_unitType.coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_unitType.corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues(2);
%                     end
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
%                 for iClass = 1:length(Class)
%                     GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).CountOverlap(iGLM,iGLM2) = sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM) == 1 & GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred2}).(Class{iClass})(:,iGLM2) == 1); 
%                     GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).Proportion(iGLM,iGLM2) = GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).CountOverlap(iGLM,iGLM2) / sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM));
%                     GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).PropTotal(iGLM,iGLM2) = GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).CountOverlap(iGLM,iGLM2) / Total_count(iClass);
%                     
%                 [GLM_unitType.corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr GLM_unitType.corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues] = corrcoef(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass})(:,iGLM),GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{iPred2}).(Class{iClass})(:,iGLM2));
%                     GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).RecodeCorr(iGLM,iGLM2) = GLM_unitType.corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorr(2);
%                     GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_unitType.corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).RecodeCorrPvalues(2);
%                 end
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
%     for iClass = 1:length(Class)
%             GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Overlap = [];
%              GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Proportion = [];
%               GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).PropTotal = [];
%               GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr = [];
%     GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue = [];
%             end
    for iPred = 1:3
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue = [];
%         for iClass = 1:length(Class)
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllOverlap = [];
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllProportion = [];
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPropTotal = [];
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredOverlap = [];
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredProportion = [];
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredPropTotal = [];
%         
%         GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Corr = [];
%         GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Pvalue = [];
%         GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Corr = [];
%         GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Pvalue = [];
%         end
        if iEpoch == 1
            GLM_coeff.summary.(Predictors{iPred}).Corr = [];
            GLM_coeff.summary.(Predictors{iPred}).Pvalue = [];
%             for iClass = 1:length(Class)
%             GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Overlap = [];
%              GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Proportion = [];
%               GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).PropTotal = [];
%               
%               GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr = [];
%             GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue = [];
%             end
        end
        for iEpoch2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(start:11,start:11));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(start:11,start:11));
%         for iClass = 1:length(Class)
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllOverlap = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllOverlap,GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).CountOverlap(start:11,start:11));
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllProportion = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllProportion,GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).Proportion(start:11,start:11));
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPropTotal = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPropTotal,GLM_unitType.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).PropTotal(start:11,start:11));
%         
%         GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Corr = cat(2,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Corr,GLM_unitType.coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).RecodeCorr(start:11,start:11));
%             GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Pvalue = cat(2,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Pvalue,GLM_unitType.coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Class{iClass}).RecodeCorrPvalues(start:11,start:11));
%         end
        end
        GLM_coeff.summary.(Predictors{iPred}).Corr = cat(1,GLM_coeff.summary.(Predictors{iPred}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr);
        GLM_coeff.summary.(Predictors{iPred}).Pvalue = cat(1,GLM_coeff.summary.(Predictors{iPred}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue);
%         for iClass = 1:length(Class)
%             GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Overlap = cat(1,GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Overlap,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllOverlap);
%              GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Proportion = cat(1,GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).Proportion,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllProportion);
%               GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).PropTotal = cat(1,GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).PropTotal,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPropTotal);
%         
%               GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr = cat(1,GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Corr);
%         GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue = cat(1,GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).ALL_Pvalue);
%         end
        for iPred2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr(start:11,start:11));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues(start:11,start:11));
%         for iClass = 1:length(Class)
%             GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredOverlap = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredOverlap,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).CountOverlap(start:11,start:11));
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredProportion = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredProportion,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).Proportion(start:11,start:11));
%         GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredPropTotal = cat(2,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredPropTotal,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).PropTotal(start:11,start:11));
%        
%          GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Corr = cat(2,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Corr,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).RecodeCorr(start:11,start:11));
%             GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Pvalue = cat(2,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Pvalue,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(Class{iClass}).RecodeCorrPvalues(start:11,start:11));       
%         end
        end
        GLM_coeff.summary.(Epoch{iEpoch}).Corr = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr);
        GLM_coeff.summary.(Epoch{iEpoch}).Pvalue = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue);
%         for iClass = 1:length(Class)
%             GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Overlap = cat(1,GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Overlap,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredOverlap);
%             GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Proportion = cat(1,GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).Proportion,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredProportion);
%             GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).PropTotal = cat(1,GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).PropTotal,GLM_unitType.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPredPropTotal);
%        
%         GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr = cat(1,GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Corr);
%         GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue = cat(1,GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue,GLM_unitType.coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Class{iClass}).AllPred_Pvalue);        
%         end
    end
end

%% overlap among all cue features or epochs
for iEpoch = 1:3
AllOverlap.(Epoch{iEpoch}) = sum(GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{1}).(Class{1})(:,6) == 1 & GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{2}).(Class{1})(:,6) == 1 & GLM_unitType.RAW.(Epoch{iEpoch}).(Predictors{3}).(Class{1})(:,6) == 1); 
end

for iPred = 1:3
AllOverlap.(Predictors{iPred}) = sum(GLM_unitType.RAW.(Epoch{1}).(Predictors{iPred}).(Class{1})(:,6) == 1 & GLM_unitType.RAW.(Epoch{2}).(Predictors{iPred}).(Class{1})(:,6) == 1 & GLM_unitType.RAW.(Epoch{3}).(Predictors{iPred}).(Class{1})(:,6) == 1); 
end

%%
% rColorMap = [linspace(233/255, 255/255, 183),linspace(255/255, 161/255, 73)];
% gColorMap = [linspace(163/255, 255/255, 183),linspace(255/255, 215/255, 73)];
% bColorMap = [linspace(201/255, 255/255, 183),linspace(255/255, 106/255, 73)];
% colorMap = [rColorMap; gColorMap; bColorMap]';
% 
% mincolor = .01;
% maxcolor = .2;
% figure
% for iPred = 1:3
%     %     subplot(2,3,iPred)
%     figure
%     GLM_coeff.summary.(Predictors{iPred}).Pvalue(GLM_coeff.summary.(Predictors{iPred}).Pvalue == 0)=NaN;
%     
%     heatmap(GLM_coeff.summary.(Predictors{iPred}).Pvalue,[],[],'%0.2f','ColorMap', flipud(colorMap),...% 'Colorbar',true, ...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0])%...
%     %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
%     % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
%     hold on
%     plot([6.5 6.5],[0.5 18.5],'k'); plot([.5 18.5],[6.5 6.5],'k');
%     plot([12.5 12.5],[0.5 18.5],'k'); plot([.5 18.5],[12.5 12.5],'k');
%     title(Predictors{iPred})
% end
% 
% %%
% for iEpoch = 1:3
%     %     subplot(2,3,iEpoch+3)
%     figure
%     GLM_coeff.summary.(Epoch{iEpoch}).Pvalue(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue == 0)=NaN;
%     
%     heatmap(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue,[],[],'%0.2f','ColorMap', flipud(colorMap),...% 'Colorbar',true, ...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0])%...
%     %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
%     % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
%     hold on
%     plot([6.5 6.5],[0.5 18.5],'k'); plot([.5 18.5],[6.5 6.5],'k');
%     plot([12.5 12.5],[0.5 18.5],'k'); plot([.5 18.5],[12.5 12.5],'k');
%     title(Epoch{iEpoch})
% end

%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};
graph_title = {'Identity coding across task epochs' 'Location coding across task epochs' 'Outcome coding across task epochs'};
colors.Modality = {[49/255 163/255 84/255] [255/255 207/255 250/255] [.8 .8 .8] ...
    [255/255 207/255 250/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
    [.8 .8 .8] [49/255 163/255 84/255] [49/255 163/255 84/255]};
 colors.Location = {[49/255 163/255 84/255] [.8 .8 .8] [49/255 163/255 84/255] ...
    [.8 .8 .8] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
    [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};
 colors.Outcome = {[49/255 163/255 84/255] [255/255 207/255 250/255] [255/255 207/255 250/255] ...
    [255/255 207/255 250/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
    [255/255 207/255 250/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};


rColorMap = [linspace(253/255, 255/255, 45),linspace(255/255, 49/255, 211)]; %77 253
gColorMap = [linspace(224/255, 255/255, 45),linspace(255/255, 163/255, 211)]; %146 224 49,163,84
bColorMap = [linspace(239/255, 255/255, 45),linspace(255/255, 84/255, 211)]; %33 239
colorMap = [rColorMap; gColorMap; bColorMap]';

labels = {'0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_cue = {'Cue identity','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue location','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_start = {'Cue-onset','','','','','','Nosepoke','','','','','','Outcome','','','','',''};
mincolor = -.2;
maxcolor = .95%.7; %.8;
% figure
for iPred = 1:3
%         subplot(2,3,iPred)
    figure
%     GLM_coeff.summary.(Predictors{iPred}).Corr(GLM_coeff.summary.(Predictors{iPred}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Predictors{iPred}).Corr(GLM_coeff.summary.(Predictors{iPred}).Pvalue > .05)=NaN;
    
    heatmap(GLM_coeff.summary.(Predictors{iPred}).Corr,[],[],[],'ColorMap',colorMap,'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [.6 .6 .6], 'TickAngle',90, 'ShowAllTicks', true)%...
    %    'RowLabels', {'onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on

   
        plot([6.27 6.27],[0.7 6.25],'color',colors.(Predictors{iPred}){1},'LineWidth',12); plot([0.7 0.7],[0.7 6.25],'color',colors.(Predictors{iPred}){1},'LineWidth',12);
    plot([0.7 6.27],[6.25 6.25],'color',colors.(Predictors{iPred}){1},'LineWidth',12); plot([0.7 6.27],[0.7 0.7],'color',colors.(Predictors{iPred}){1},'LineWidth',12);
    plot([12.27 12.27],[0.7 6.25],'color',colors.(Predictors{iPred}){2},'LineWidth',12); plot([6.73 6.73],[0.7 6.25],'color',colors.(Predictors{iPred}){2},'LineWidth',12);
    plot([6.73 12.27],[6.25 6.25],'color',colors.(Predictors{iPred}){2},'LineWidth',12); plot([6.73 12.27],[0.7 0.7],'color',colors.(Predictors{iPred}){2},'LineWidth',12);
  plot([18.3 18.3],[0.7 6.25],'color',colors.(Predictors{iPred}){3},'LineWidth',12); plot([12.73 12.73],[0.7 6.25],'color',colors.(Predictors{iPred}){3},'LineWidth',12);
    plot([12.73 18.3],[6.25 6.25],'color',colors.(Predictors{iPred}){3},'LineWidth',12); plot([12.73 18.3],[0.7 0.7],'color',colors.(Predictors{iPred}){3},'LineWidth',12);
    
     plot([6.27 6.27],[6.75 12.25],'color',colors.(Predictors{iPred}){4},'LineWidth',12); plot([0.7 0.7],[6.75 12.25],'color',colors.(Predictors{iPred}){4},'LineWidth',12);
    plot([0.7 6.27],[12.25 12.25],'color',colors.(Predictors{iPred}){4},'LineWidth',12); plot([0.7 6.27],[6.75 6.75],'color',colors.(Predictors{iPred}){4},'LineWidth',12);
    plot([12.27 12.27],[6.75 12.25],'color',colors.(Predictors{iPred}){5},'LineWidth',12); plot([6.73 6.73],[6.75 12.25],'color',colors.(Predictors{iPred}){5},'LineWidth',12);
    plot([6.73 12.27],[12.25 12.25],'color',colors.(Predictors{iPred}){5},'LineWidth',12); plot([6.73 12.27],[6.75 6.75],'color',colors.(Predictors{iPred}){5},'LineWidth',12);
  plot([18.3 18.3],[6.75 12.25],'color',colors.(Predictors{iPred}){6},'LineWidth',12); plot([12.73 12.73],[6.75 12.25],'color',colors.(Predictors{iPred}){6},'LineWidth',12);
    plot([12.73 18.3],[12.25 12.25],'color',colors.(Predictors{iPred}){6},'LineWidth',12); plot([12.73 18.3],[6.75 6.75],'color',colors.(Predictors{iPred}){6},'LineWidth',12);
    
         plot([6.27 6.27],[12.75 18.3],'color',colors.(Predictors{iPred}){7},'LineWidth',12); plot([0.7 0.7],[12.75 18.3],'color',colors.(Predictors{iPred}){7},'LineWidth',12);
    plot([0.7 6.27],[18.3 18.3],'color',colors.(Predictors{iPred}){7},'LineWidth',12); plot([0.7 6.27],[12.75 12.75],'color',colors.(Predictors{iPred}){7},'LineWidth',12);
    plot([12.27 12.27],[12.75 18.3],'color',colors.(Predictors{iPred}){8},'LineWidth',12); plot([6.73 6.73],[12.75 18.3],'color',colors.(Predictors{iPred}){8},'LineWidth',12);
    plot([6.73 12.27],[18.3 18.3],'color',colors.(Predictors{iPred}){8},'LineWidth',12); plot([6.73 12.27],[12.75 12.75],'color',colors.(Predictors{iPred}){8},'LineWidth',12);
  plot([18.3 18.3],[12.75 18.3],'color',colors.(Predictors{iPred}){9},'LineWidth',12); plot([12.73 12.73],[12.75 18.3],'color',colors.(Predictors{iPred}){9},'LineWidth',12);
    plot([12.73 18.3],[18.3 18.3],'color',colors.(Predictors{iPred}){9},'LineWidth',12); plot([12.73 18.3],[12.75 12.75],'color',colors.(Predictors{iPred}){9},'LineWidth',12);

             plot([6.5 6.5],[0.56 19.49],'k','LineWidth',4); plot([-0.51 18.44],[6.5 6.5],'k','LineWidth',4);
    plot([12.5 12.5],[0.56 19.49],'k','LineWidth',4); plot([-0.51 18.44],[12.5 12.5],'k','LineWidth',4);
    % title(Predictors{iPred})
    % xlabel('cue onset')
     xlim([0.5 18.5])
     ylim([0.5 18.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
     y = ylabel('Outcome                   Nosepoke                   Cue-onset');
     x = xlabel('Cue-onset                        Nosepoke                         Outcome');
     ax = gca;
     ax.Clipping = 'off';
     title(graph_title{iPred})
end

%%
graph_title = {'Coding of cue features at cue-onset' 'Coding of cue features at nosepoke' 'Coding of cue features at outcome receipt'};
colors.cueon = {[49/255 163/255 84/255] [.8 .8 .8] [255/255 207/255 250/255] ...
    [.8 .8 .8] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
    [255/255 207/255 250/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};
 colors.NP = {[49/255 163/255 84/255] [.8 .8 .8] [255/255 207/255 250/255] ...
    [.8 .8 .8] [49/255 163/255 84/255] [.8 .8 .8] ...
   [255/255 207/255 250/255] [.8 .8 .8] [49/255 163/255 84/255]};
 colors.outcome = {[49/255 163/255 84/255] [.8 .8 .8] [255/255 207/255 250/255] ...
    [.8 .8 .8] [49/255 163/255 84/255] [255/255 207/255 250/255] ...
    [255/255 207/255 250/255] [255/255 207/255 250/255] [49/255 163/255 84/255]};

for iEpoch = 1:3
%         subplot(2,3,iEpoch+3)
    figure
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue > .05)=NaN;
    
    heatmap(GLM_coeff.summary.(Epoch{iEpoch}).Corr,[],[],[],'ColorMap',colorMap,'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0.6 0.6 0.6], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
   
       
    
        plot([6.27 6.27],[0.7 6.25],'color',colors.(Epoch{iEpoch}){1},'LineWidth',12); plot([0.7 0.7],[0.7 6.25],'color',colors.(Epoch{iEpoch}){1},'LineWidth',12);
    plot([0.7 6.27],[6.25 6.25],'color',colors.(Epoch{iEpoch}){1},'LineWidth',12); plot([0.7 6.27],[0.7 0.7],'color',colors.(Epoch{iEpoch}){1},'LineWidth',12);
    plot([12.27 12.27],[0.7 6.25],'color',colors.(Epoch{iEpoch}){2},'LineWidth',12); plot([6.73 6.73],[0.7 6.25],'color',colors.(Epoch{iEpoch}){2},'LineWidth',12);
    plot([6.73 12.27],[6.25 6.25],'color',colors.(Epoch{iEpoch}){2},'LineWidth',12); plot([6.73 12.27],[0.7 0.7],'color',colors.(Epoch{iEpoch}){2},'LineWidth',12);
  plot([18.3 18.3],[0.7 6.25],'color',colors.(Epoch{iEpoch}){3},'LineWidth',12); plot([12.73 12.73],[0.7 6.25],'color',colors.(Epoch{iEpoch}){3},'LineWidth',12);
    plot([12.73 18.3],[6.25 6.25],'color',colors.(Epoch{iEpoch}){3},'LineWidth',12); plot([12.73 18.3],[0.7 0.7],'color',colors.(Epoch{iEpoch}){3},'LineWidth',12);
    
     plot([6.27 6.27],[6.75 12.25],'color',colors.(Epoch{iEpoch}){4},'LineWidth',12); plot([0.7 0.7],[6.75 12.25],'color',colors.(Epoch{iEpoch}){4},'LineWidth',12);
    plot([0.7 6.27],[12.25 12.25],'color',colors.(Epoch{iEpoch}){4},'LineWidth',12); plot([0.7 6.27],[6.75 6.75],'color',colors.(Epoch{iEpoch}){4},'LineWidth',12);
    plot([12.27 12.27],[6.75 12.25],'color',colors.(Epoch{iEpoch}){5},'LineWidth',12); plot([6.73 6.73],[6.75 12.25],'color',colors.(Epoch{iEpoch}){5},'LineWidth',12);
    plot([6.73 12.27],[12.25 12.25],'color',colors.(Epoch{iEpoch}){5},'LineWidth',12); plot([6.73 12.27],[6.75 6.75],'color',colors.(Epoch{iEpoch}){5},'LineWidth',12);
  plot([18.3 18.3],[6.75 12.25],'color',colors.(Epoch{iEpoch}){6},'LineWidth',12); plot([12.73 12.73],[6.75 12.25],'color',colors.(Epoch{iEpoch}){6},'LineWidth',12);
    plot([12.73 18.3],[12.25 12.25],'color',colors.(Epoch{iEpoch}){6},'LineWidth',12); plot([12.73 18.3],[6.75 6.75],'color',colors.(Epoch{iEpoch}){6},'LineWidth',12);
    
         plot([6.27 6.27],[12.75 18.3],'color',colors.(Epoch{iEpoch}){7},'LineWidth',12); plot([0.7 0.7],[12.75 18.3],'color',colors.(Epoch{iEpoch}){7},'LineWidth',12);
    plot([0.7 6.27],[18.3 18.3],'color',colors.(Epoch{iEpoch}){7},'LineWidth',12); plot([0.7 6.27],[12.75 12.75],'color',colors.(Epoch{iEpoch}){7},'LineWidth',12);
    plot([12.27 12.27],[12.75 18.3],'color',colors.(Epoch{iEpoch}){8},'LineWidth',12); plot([6.73 6.73],[12.75 18.3],'color',colors.(Epoch{iEpoch}){8},'LineWidth',12);
    plot([6.73 12.27],[18.3 18.3],'color',colors.(Epoch{iEpoch}){8},'LineWidth',12); plot([6.73 12.27],[12.75 12.75],'color',colors.(Epoch{iEpoch}){8},'LineWidth',12);
  plot([18.3 18.3],[12.75 18.3],'color',colors.(Epoch{iEpoch}){9},'LineWidth',12); plot([12.73 12.73],[12.75 18.3],'color',colors.(Epoch{iEpoch}){9},'LineWidth',12);
    plot([12.73 18.3],[18.3 18.3],'color',colors.(Epoch{iEpoch}){9},'LineWidth',12); plot([12.73 18.3],[12.75 12.75],'color',colors.(Epoch{iEpoch}){9},'LineWidth',12);

    plot([6.5 6.5],[0.56 19.49],'k','LineWidth',4); plot([-.51 18.44],[6.5 6.5],'k','LineWidth',4);
    plot([12.5 12.5],[0.56 19.49],'k','LineWidth',4); plot([-.51 18.44],[12.5 12.5],'k','LineWidth',4);
   
   % title(Epoch{iEpoch})
        xlim([0.5 18.5])
     ylim([0.5 18.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
%      y = ylabel('Cue outcome              Cue location               Cue identity');
%      x = xlabel('Cue identity                      Cue location                     Cue outcome');
y = ylabel('Outcome coding          Location coding           Identity coding');
     x = xlabel('Identity coding               Location coding              Outcome coding');
     ax = gca;
     ax.Clipping = 'off';
     title(graph_title{iEpoch})
end

%%
Fade = [1 .6 .475 .35 .225 .1];
matrix_start = [1 7 13];
joint_input = [.1 .2 .3 .4 .5];
separate_input = [-.2 -1.55 -.1 -.075 -.03];
ind_input = [0 -.04 .05 -.025 .035];
schematic_data(1:18,1:18) = NaN; 
for iFade = 1:6
    schematic_data(iFade:18+1:end) = Fade(iFade);
    schematic_data(21-iFade:18+1:end) = Fade(iFade);
end
schematic_data_joint = schematic_data;
schematic_data_sep = schematic_data;
schematic_data_ind = schematic_data;
for iBlockRow = 1:3
    for iBlockCol = 1:3
        if iBlockRow ~= iBlockCol
joint_data = datasample(joint_input,36);
% joint_data = .1:.1:3.6;
joint_data2 = cat(1,joint_data(1:6),joint_data(7:12),joint_data(13:18),joint_data(19:24),joint_data(25:30),joint_data(31:36));
separate_data = datasample(separate_input,36);
NaNs_data(1:6,1:6) = NaN;
% separate_data = .1:.1:3.6;
separate_data2 = cat(1,separate_data(1:6),separate_data(7:12),separate_data(13:18),separate_data(19:24),separate_data(25:30),separate_data(31:36));
ind_data = datasample(ind_input,36);
% ind_data = .1:.1:3.6;
ind_data2 = cat(1,ind_data(1:6),ind_data(7:12),ind_data(13:18),ind_data(19:24),ind_data(25:30),ind_data(31:36));

schematic_data_joint(matrix_start(iBlockRow):matrix_start(iBlockRow)+5,matrix_start(iBlockCol):matrix_start(iBlockCol)+5) = joint_data2;
schematic_data_sep(matrix_start(iBlockRow):matrix_start(iBlockRow)+5,matrix_start(iBlockCol):matrix_start(iBlockCol)+5) = separate_data2;
schematic_data_ind(matrix_start(iBlockRow):matrix_start(iBlockRow)+5,matrix_start(iBlockCol):matrix_start(iBlockCol)+5) = ind_data2;
  schematic_data(matrix_start(iBlockRow):matrix_start(iBlockRow)+5,matrix_start(iBlockCol):matrix_start(iBlockCol)+5) = NaNs_data;
        end
        end
end
%%
for iPlot = 1:3
    switch iPlot
        case 1
            plot_schematic = schematic_data_ind;
        case 2
            plot_schematic = schematic_data_joint;
        case 3
            plot_schematic = schematic_data_sep;
    end

figure
    heatmap(plot_schematic,[],[],[],'ColorMap',colorMap,...%'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0.6 0.6 0.6], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.56 18.49],'k','LineWidth',2); plot([.51 18.44],[6.5 6.5],'k','LineWidth',2);
    plot([12.5 12.5],[0.56 18.49],'k','LineWidth',2); plot([.51 18.44],[12.5 12.5],'k','LineWidth',2);
    % title(Epoch{iEpoch})
        xlim([0.5 18.5])
     ylim([0.5 18.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
%      y = ylabel('Cue outcome              Cue location               Cue identity');
%      x = xlabel('Cue identity                      Cue location                     Cue outcome');
% y = ylabel('Outcome coding          Location coding           Identity coding');
%      x = xlabel('Identity coding               Location coding              Outcome coding');
     ax = gca;
     ax.Clipping = 'off';
%      title(graph_title{iEpoch})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cue type other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};

rColorMap = [linspace(255/255, 161/255, 256)]; %77
gColorMap = [linspace(255/255, 215/255, 256)]; %146
bColorMap = [linspace(255/255, 106/255, 256)]; %33
colorMap = [rColorMap; gColorMap; bColorMap]';

labels = {'Cue-onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Nosepoke','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_cue = {'Cue identity','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue location','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
% mincolor = -.4;
% maxcolor = .7; %.8;
Type = {'Overlap' 'Proportion' 'PropTotal'};
for iType = 1%1:3
% figure
iPlot = 1;
for iPred = 1:3
    for iClass = 1:3
%          subplot(3,3,iPlot)
         iPlot = iPlot+1;
    figure
%     GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue == 0)=NaN;
%  GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr(eye(size(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue)) == 1)=NaN;
%     GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).(Type{iType})(GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).(Type{iType}) < .4)=NaN;
    
    heatmap(GLM_unitType.summary.(Predictors{iPred}).(Class{iClass}).(Type{iType}),labels,labels,'%0.0f','ColorMap',colorMap,'NaNColor', [0 0 0])%,... 'Colorbar',true, ...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[6.5 6.5],'k','LineWidth',2);
    plot([12.5 12.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[12.5 12.5],'k','LineWidth',2);
    title(strcat(Predictors{iPred},' (',Class{iClass},')'))
    % xlabel('cue onset')
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
    end
end
end

%%
for iType = 1%1:3
% figure
iPlot = 1;
for iEpoch = 1:3
    for iClass = 1%:3
%         subplot(3,3,iPlot)
    figure
iPlot = iPlot + 1;
%     GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue < .0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)=NaN;
% GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(eye(size(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue)) == 1)=NaN;
%     GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue > .05)=NaN;
% GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).(Type{iType})(GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).(Type{iType}) < .4)=NaN;    

        heatmap(GLM_unitType.summary.(Epoch{iEpoch}).(Class{iClass}).(Type{iType}),labels_cue,labels_cue,'%0.0f','ColorMap',colorMap,'NaNColor', [0 0 0])%... 'Colorbar',true, ...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %     heatmap(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr,labels_cue,labels_cue,'%0.1f','ColorMap',colorMap,... 'Colorbar',true, ...

    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[6.5 6.5],'k','LineWidth',2);
    plot([12.5 12.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[12.5 12.5],'k','LineWidth',2);
    title(strcat(Epoch{iEpoch},' (',Class{iClass},')'))
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
    end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CELL TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};

rColorMap = [linspace(233/255, 255/255, 73),linspace(255/255, 161/255, 183)]; %77
gColorMap = [linspace(163/255, 255/255, 73),linspace(255/255, 215/255, 183)]; %146
bColorMap = [linspace(201/255, 255/255, 73),linspace(255/255, 106/255, 183)]; %33
colorMap = [rColorMap; gColorMap; bColorMap]';

labels = {'Cue-onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Nosepoke','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_cue = {'Cue identity','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue location','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
mincolor = -.4;
maxcolor = .7; %.8;
figure
iPlot = 1;
for iPred = 1:3
    for iClass = 1:3
         subplot(3,3,iPlot)
         iPlot = iPlot+1;
%     figure
%     GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue == 0)=NaN;
 GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr(eye(size(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue)) == 1)=NaN;
    GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Pvalue > .05)=NaN;
    
    heatmap(GLM_unitType.coeff.summary.(Predictors{iPred}).(Class{iClass}).Corr,[],[],[],'ColorMap',colorMap,... 'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[6.5 6.5],'k','LineWidth',2);
    plot([12.5 12.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[12.5 12.5],'k','LineWidth',2);
    title(strcat(Predictors{iPred},' (',Class{iClass},')'))
    % xlabel('cue onset')
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
    end
end

%%
figure
iPlot = 1;
for iEpoch = 1:3
    for iClass = 1:3
        subplot(3,3,iPlot)
%     figure
iPlot = iPlot + 1;
%     GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue < .0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)=NaN;
GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(eye(size(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue)) == 1)=NaN;
    GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Pvalue > .05)=NaN;
    
        heatmap(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr,[],[],[],'ColorMap',colorMap,... 'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %     heatmap(GLM_unitType.coeff.summary.(Epoch{iEpoch}).(Class{iClass}).Corr,labels_cue,labels_cue,'%0.1f','ColorMap',colorMap,... 'Colorbar',true, ...

    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[6.5 6.5],'k','LineWidth',2);
    plot([12.5 12.5],[0.51 18.49],'k','LineWidth',2); plot([.51 18.49],[12.5 12.5],'k','LineWidth',2);
    title(strcat(Epoch{iEpoch},' (',Class{iClass},')'))
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
    end
end