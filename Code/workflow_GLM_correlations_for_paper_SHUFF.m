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
                for iShuff = 1:num_Shuffs
                for iPred = 1:length(Predictors)                 
                    count = 1;
                    count_pres = 1;
                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(1:133,iGLM) = 0;
                    for iCMod = 1:length(cue_mod)
                        if iCMod <= length(ALL_matrix.(cat(2,'shuff_',num2str(iShuff))))
                            if cue_mod(iCMod) == 1
                                if cue_presence(count_pres) == 1
                                if ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(iCMod,iPred) == 1
                                    GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(count,iGLM) = 1;
                                end
                                count = count + 1;
                                end
                                count_pres = count_pres + 1;
                            end
                        end
                    end
                end
                end
            end
        end
    end
end

%% across epochs
rng('shuffle')
num_Shuffs = 100;
num_GLM = 11; %6 or 11 depending if want all time windows
% Class = {'ALL' 'MSN' 'FSI' 'MSN_Inc' 'MSN_Dec' 'FSI_Inc' 'FSI_Dec'};
% Total_count = [sum(cue_mod) sum(cue_MSN) sum(cue_FSI) sum(cue_MSN_Inc) sum(cue_MSN_Dec) sum(cue_FSI_Inc) sum(cue_FSI_Dec)];

for iEpoch = 1:3
    for iEpoch2 = 1:3
        for iGLM = 1:num_GLM
            for iGLM2 = 1:num_GLM
                disp(cat(2,'Epochs ',num2str(iEpoch),' (',num2str(iGLM),') & ',num2str(iEpoch2),' (',num2str(iGLM2),')'))
                for iPred = 1:3
%                     Shuff_comp = randperm(100);
                    for iShuff = 1:num_Shuffs
%                         iShuff2 = Shuff_comp(iShuff);
%                         Epoch1_shuff = [];
%                         Epoch2_shuff = [];
%                         Epoch1_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM)),'Replace',false);
%                         Epoch2_shuff = datasample(GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2),length(GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred})(:,iGLM2)),'Replace',false);                          
%                     [GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(Epoch1_shuff,Epoch2_shuff);
                    [GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(:,iGLM),GLM_recode.(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(:,iGLM2));
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr(2);
                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(2);
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
%                     Shuff_comp = randperm(100);
                    for iShuff = 1:num_Shuffs
%                         iShuff2 = Shuff_comp(iShuff)
%                         Pred1_shuff = [];
%                         Pred2_shuff = [];
%                         Pred1_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred})(:,iGLM)),'Replace',false);
%                         Pred2_shuff = datasample(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2})(:,iGLM2),length(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2})(:,iGLM2)),'Replace',false);   
%                     [GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(Pred1_shuff,Pred2_shuff);
                    [GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues] = corrcoef(GLM_recode.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff)))(:,iGLM),GLM_recode.(Epoch{iEpoch}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff)))(:,iGLM2));
                    GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorr(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorr(2);
                    GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(iGLM,iGLM2) = GLM_corr.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(strcat('t',num2str(iGLM))).(strcat('t',num2str(iGLM2))).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(2);
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
    
    GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr = [];
    GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue = [];
    for iPred = 1:3
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr = [];
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue = [];
        if iEpoch == 1
            GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr = [];
            GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue = [];
        end
        for iEpoch2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorr(start:num_GLM,start:num_GLM));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(start:num_GLM,start:num_GLM));
        end
        GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr = cat(1,GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Corr);
        GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue = cat(1,GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).ALL_Pvalue);
        for iPred2 = 1:3
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorr(start:num_GLM,start:num_GLM));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).(cat(2,'s',num2str(iShuff))).RecodeCorrPvalues(start:num_GLM,start:num_GLM));
        end
        GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Corr);
        GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).AllPred_Pvalue);
    end
end
end

GLM_coeff_SHUFF = GLM_coeff;
%% Z-scores away from shuffle (value - shuff mean / shuff std)
%%%% for MEAN of Corr %%%%
matrix_start = [1 7 13];
num_Shuffs = 100;
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
    %%

%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};
graph_title = {'Identity coding across task epochs' 'Location coding across task epochs' 'Outcome coding across task epochs'};

rColorMap = [linspace(253/255, 255/255, 45),linspace(255/255, 49/255, 211)]; %77 253
gColorMap = [linspace(224/255, 255/255, 45),linspace(255/255, 163/255, 211)]; %146 224 49,163,84
bColorMap = [linspace(239/255, 255/255, 45),linspace(255/255, 84/255, 211)]; %33 239
colorMap = [rColorMap; gColorMap; bColorMap]';

labels = {'0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_cue = {'Cue identity','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue location','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_start = {'Cue-onset','','','','','','Nosepoke','','','','','','Outcome','','','','',''};
mincolor = -.2;
maxcolor = .95%.7; %.8;
figure
% for iShuff = 1:num_Shuffs
for iPred = 1:3
        subplot(2,3,iPred)
%     figure
%     GLM_coeff.summary.(Predictors{iPred}).Corr(GLM_coeff.summary.(Predictors{iPred}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr(GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Pvalue > .05)=NaN;
     heatmap(Table.(Predictors{iPred}).Zscore_recode,[],[],[],'ColorMap',colorMap,... %'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [.6 .6 .6], 'TickAngle',90, 'ShowAllTicks', true)%...
   
%     heatmap(GLM_coeff.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr,[],[],[],'ColorMap',colorMap,'Colorbar',true, ...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [.6 .6 .6], 'TickAngle',90, 'ShowAllTicks', true)%...
%     %    'RowLabels', {'onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([2.5 2.5],[0.5 3.5],'k','LineWidth',1); plot([0.5 3.5],[2.5 2.5],'k','LineWidth',1);
    plot([1.5 1.5],[0.5 3.5],'k','LineWidth',1); plot([0.5 3.5],[1.5 1.5],'k','LineWidth',1);
    % title(Predictors{iPred})
    % xlabel('cue onset')
     xlim([0.5 3.5])
     ylim([0.5 3.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
%      y = ylabel('Outcome                   Nosepoke                   Cue-onset');
%      x = xlabel('Cue-onset                        Nosepoke                         Outcome');
     ax = gca;
     ax.Clipping = 'off';
     title(graph_title{iPred})
end
end

%%
graph_title = {'Cue features at cue-onset' 'Cue features at nosepoke' 'Cue features at outcome receipt'};


for iEpoch = 1:3
        subplot(2,3,iEpoch+3)
%     figure
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue > .05)=NaN;
    
    heatmap(Table.(Epoch{iEpoch}).Zscore_recode,[],[],[],'ColorMap',colorMap,... %'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0.6 0.6 0.6], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
 plot([2.5 2.5],[0.5 3.5],'k','LineWidth',1); plot([0.5 3.5],[2.5 2.5],'k','LineWidth',1);
    plot([1.5 1.5],[0.5 3.5],'k','LineWidth',1); plot([0.5 3.5],[1.5 1.5],'k','LineWidth',1);
    % title(Predictors{iPred})
    % xlabel('cue onset')
     xlim([0.5 3.5])
     ylim([0.5 3.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
%      y = ylabel('Cue outcome              Cue location               Cue identity');
%      x = xlabel('Cue identity                      Cue location                     Cue outcome');
% y = ylabel('Outcome coding          Location coding           Identity coding');
%      x = xlabel('Identity coding               Location coding              Outcome coding');
     ax = gca;
     ax.Clipping = 'off';
     title(graph_title{iEpoch})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Z-scores away from shuffle (value - shuff mean / shuff std)
%%%% for individual Corrs %%%%
% matrix_start = [1 7 13];
for iShuff = 1:num_Shuffs
    for iEpoch = 1:3
        for iRow = 1:18
            for iCol = 1:18
                Table_ind.Shuffs.(Epoch{iEpoch}){iRow,iCol}(iShuff) = GLM_coeff_SHUFF.summary.(Epoch{iEpoch}).(cat(2,'s',num2str(iShuff))).Corr(iRow,iCol);
            end
        end
    end
    
    for iPred = 1:3
        for iRow = 1:18
            for iCol = 1:18
                Table_ind.Shuffs.(Predictors{iPred}){iRow,iCol}(iShuff) = GLM_coeff_SHUFF.summary.(Predictors{iPred}).(cat(2,'s',num2str(iShuff))).Corr(iRow,iCol);
            end
        end
    end
end

for iEpoch = 1:3
    for iRow = 1:18
        for iCol = 1:18
            Table_ind.(Epoch{iEpoch}).Data(iRow,iCol) = GLM_coeff.summary.(Epoch{iEpoch}).Corr(iRow,iCol);
            Table_ind.(Epoch{iEpoch}).ShuffMEAN(iRow,iCol) = mean(Table_ind.Shuffs.(Epoch{iEpoch}){iRow,iCol});
            Table_ind.(Epoch{iEpoch}).ShuffSTD(iRow,iCol) = std(Table_ind.Shuffs.(Epoch{iEpoch}){iRow,iCol});
        Table_ind.(Epoch{iEpoch}).Zscore(iRow,iCol) = (Table_ind.(Epoch{iEpoch}).Data(iRow,iCol) - Table_ind.(Epoch{iEpoch}).ShuffMEAN(iRow,iCol)) / Table_ind.(Epoch{iEpoch}).ShuffSTD(iRow,iCol);
       if Table_ind.(Epoch{iEpoch}).Zscore(iRow,iCol) > 1.96
                Table_ind.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = 1;
       elseif Table_ind.(Epoch{iEpoch}).Zscore(iRow,iCol) < -1.96
                Table_ind.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = -1;
            else
                Table_ind.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = 0;
       end
             if iRow == iCol
                 Table_ind.(Epoch{iEpoch}).Zscore_recode(iRow,iCol) = 0;
            end
        end
    end
    Table_ind.MAX(iEpoch,:) = max(Table_ind.(Epoch{iEpoch}).Zscore);
end

for iPred = 1:3
    for iRow = 1:18
        for iCol = 1:18
            Table_ind.(Predictors{iPred}).Data(iRow,iCol) = GLM_coeff.summary.(Predictors{iPred}).Corr(iRow,iCol);
            Table_ind.(Predictors{iPred}).ShuffMEAN(iRow,iCol) = mean(Table_ind.Shuffs.(Predictors{iPred}){iRow,iCol});
            Table_ind.(Predictors{iPred}).ShuffSTD(iRow,iCol) = std(Table_ind.Shuffs.(Predictors{iPred}){iRow,iCol});
            Table_ind.(Predictors{iPred}).Zscore(iRow,iCol) = (Table_ind.(Predictors{iPred}).Data(iRow,iCol) - Table_ind.(Predictors{iPred}).ShuffMEAN(iRow,iCol)) / Table_ind.(Predictors{iPred}).ShuffSTD(iRow,iCol);
       if Table_ind.(Predictors{iPred}).Zscore(iRow,iCol) > 1.96
                Table_ind.(Predictors{iPred}).Zscore_recode(iRow,iCol) = 1;
       elseif Table_ind.(Predictors{iPred}).Zscore(iRow,iCol) < -1.96
           Table_ind.(Predictors{iPred}).Zscore_recode(iRow,iCol) = -1;
       else
                Table_ind.(Predictors{iPred}).Zscore_recode(iRow,iCol) = 0;
       end
            if iRow == iCol
                 Table_ind.(Predictors{iPred}).Zscore_recode(iRow,iCol) = 0;
            end
        end
    end
    Table_ind.MAX(3+iPred,:) = max(Table_ind.(Predictors{iPred}).Zscore);
end

%% Individual corr z score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};
graph_title = {'Identity coding across task epochs' 'Location coding across task epochs' 'Outcome coding across task epochs'};

rColorMap = [linspace(253/255, 255/255, 128),linspace(255/255, 49/255, 128)]; %77 253
gColorMap = [linspace(224/255, 255/255, 128),linspace(255/255, 163/255, 128)]; %146 224 49,163,84
bColorMap = [linspace(239/255, 255/255, 128),linspace(255/255, 84/255, 128)]; %33 239
colorMap = [rColorMap; gColorMap; bColorMap]';

labels = {'0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','0.0 s','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_cue = {'Cue identity','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue location','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s','Cue outcome','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'};
labels_start = {'Cue-onset','','','','','','Nosepoke','','','','','','Outcome','','','','',''};
mincolor = -1.95;
maxcolor = 1.95%.7; %.8;
figure
for iPred = 1:3
        subplot(2,3,iPred)
%     figure
%     GLM_coeff.summary.(Predictors{iPred}).Corr(GLM_coeff.summary.(Predictors{iPred}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Predictors{iPred}).Corr(GLM_coeff.summary.(Predictors{iPred}).Pvalue > .05)=NaN;
    
    heatmap(Table_ind.(Predictors{iPred}).Zscore_recode,[],[],[],'ColorMap',colorMap,...%'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [.6 .6 .6], 'TickAngle',90, 'ShowAllTicks', true)%...
    %    'RowLabels', {'onset','+ 0.1 s','+ 0.2 s','+ 0.3 s','+ 0.4 s','+ 0.5 s'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.56 19.49],'k','LineWidth',4); plot([-0.51 18.44],[6.5 6.5],'k','LineWidth',4);
    plot([12.5 12.5],[0.56 19.49],'k','LineWidth',4); plot([-0.51 18.44],[12.5 12.5],'k','LineWidth',4);
    % title(Predictors{iPred})
    % xlabel('cue onset')
     xlim([0.5 18.5])
     ylim([0.5 18.5])
    set(gca,'FontSize',18)
    set(gcf,'Position', [10, 10, 1150, 950])
%      y = ylabel('Outcome                   Nosepoke                   Cue-onset');
%      x = xlabel('Cue-onset                        Nosepoke                         Outcome');
     ax = gca;
     ax.Clipping = 'off';
     title(graph_title{iPred})
end

%%
graph_title = {'Coding of cue features at cue-onset' 'Coding of cue features at nosepoke' 'Coding of cue features at outcome receipt'};


for iEpoch = 1:3
        subplot(2,3,iEpoch+3)
%     figure
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue == 0)=NaN;
%     GLM_coeff.summary.(Epoch{iEpoch}).Corr(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue > .05)=NaN;
    
    heatmap(Table_ind.(Epoch{iEpoch}).Zscore_recode,[],[],[],'ColorMap',colorMap,...%'Colorbar',true, ...
        'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0.6 0.6 0.6], 'TickAngle', 45, 'ShowAllTicks', true)%...
    %    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
    % set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
    hold on
    plot([6.5 6.5],[0.56 19.49],'k','LineWidth',4); plot([-.51 18.44],[6.5 6.5],'k','LineWidth',4);
    plot([12.5 12.5],[0.56 19.49],'k','LineWidth',4); plot([-.51 18.44],[12.5 12.5],'k','LineWidth',4);
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
     title(graph_title{iEpoch})
end

%% visualize histogram of shuffled data
for iEpoch = 1:3
    figure
    
    for iPlot = 1:prod(size(Table_ind.Shuffs.(Epoch{iEpoch})))
        subtightplot(18,18,iPlot)
        histfit(Table_ind.Shuffs.(Epoch{iEpoch}){iPlot},10)
        box off
        hold on
        plot([0 0],[0 25],'g')
           xlim([-.5 .5])
     ylim([0 25])
      set(gca,'XTick',[],'YTick',[]);
    end
   subtightplot(18,18,9)
   title(Epoch{iEpoch})
end

for iPred = 1:3
    figure
    
    for iPlot = 1:prod(size(Table_ind.Shuffs.(Predictors{iPred})))
        subtightplot(18,18,iPlot)
        histfit(Table_ind.Shuffs.(Predictors{iPred}){iPlot},10)
        box off
        hold on
        plot([0 0],[0 25],'g')
           xlim([-.5 .5])
     ylim([0 25])
      set(gca,'XTick',[],'YTick',[]);
    end
   subtightplot(18,18,9)
   title(Predictors{iPred})
end