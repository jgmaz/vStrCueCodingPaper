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
            if RANK.two.Trial > 975 || RANK.two.Trial < 26
                if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                    cue_mod(kk) = 1 ;
                end
            end
    end
end

%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Predictors = {'Modality' 'Location' 'Outcome'};    
for iEpoch = 1:3%:length(Epoch)
    switch iEpoch
        case 1
            % mat_files = dir('2018-03-12-GLM_cueon*');
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous' ...
            %     'ModxLoc' 'ModxOut' 'LocxOut' 'OutxApp' 'ModxLocxOut'};
            % file_order = {'2018-05-24-GLM_cueon_-0.5-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.4-recodePrevTrial.mat' '2018-05-24_cueon_-0.3-recodePrevTrial.mat' ...
            %     '2018-05-24-GLM_cueon_-0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.1-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.1-recodePrevTrial.mat' ...
            %     '2018-05-24-GLM_cueon_0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.3-recodePrevTrial.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
            file_order = {'2018-03-12-GLM_cueon_-0.5.mat' '2018-03-12-GLM_cueon_-0.4.mat' '2018-03-12-GLM_cueon_-0.3.mat' ...
                '2018-03-12-GLM_cueon_-0.2.mat' '2018-03-12-GLM_cueon_-0.1.mat' '2018-03-12-GLM_cueon_0.mat' '2018-03-12-GLM_cueon_0.1.mat' ...
                '2018-03-12-GLM_cueon_0.2.mat' '2018-03-12-GLM_cueon_0.3.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
        case 2
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-12-GLM_NP_-0.5-500bin.mat' '2018-03-12-GLM_NP_-0.4-500bin.mat' '2018-03-12-GLM_NP_-0.3-500bin.mat' ...
                '2018-03-12-GLM_NP_-0.2-500bin.mat' '2018-03-12-GLM_NP_-0.1-500bin.mat' '2018-03-12-GLM_NP_0-500bin.mat' '2018-03-12-GLM_NP_0.1-500bin.mat' ...
                '2018-03-12-GLM_NP_0.2-500bin.mat' '2018-03-12-GLM_NP_0.3-500bin.mat' '2018-03-12-GLM_NP_0.4-500bin.mat' '2018-03-12-GLM_NP_0.5-500bin.mat'};
        case 3
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-23-GLM_outcome_-0.5.mat' '2018-03-23-GLM_outcome_-0.4.mat' '2018-03-23-GLM_outcome_-0.3.mat' ...
                '2018-03-23-GLM_outcome_-0.2.mat' '2018-03-23-GLM_outcome_-0.1.mat' '2018-03-23-GLM_outcome_0.mat' '2018-03-23-GLM_outcome_0.1.mat' ...
                '2018-03-23-GLM_outcome_0.2.mat' '2018-03-23-GLM_outcome_0.3.mat' '2018-03-23-GLM_outcome_0.4.mat' '2018-03-23-GLM_outcome_0.5.mat'};
        case 4
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-12-GLM_cueoff*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-12-GLM_cueoff_-0.5.mat' '2018-03-12-GLM_cueoff_-0.4.mat' '2018-03-12-GLM_cueoff_-0.3.mat' ...
                '2018-03-12-GLM_cueoff_-0.2.mat' '2018-03-12-GLM_cueoff_-0.1.mat' '2018-03-12-GLM_cueoff_0.mat' '2018-03-12-GLM_cueoff_0.1.mat' ...
                '2018-03-12-GLM_cueoff_0.2.mat' '2018-03-12-GLM_cueoff_0.3.mat' '2018-03-12-GLM_cueoff_0.4.mat' '2018-03-12-GLM_cueoff_0.5.mat'};
    end
    
    for iGLM = 1:length(file_order)
        for iFind = 1:length(mat_files)
            if strcmp(mat_files(iFind).name,file_order{iGLM}) == 1
                load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files(iFind).name),'ALL_matrix');                                            
                           
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
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(6:11,6:11));
            GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(6:11,6:11));
        end
        GLM_coeff.summary.(Predictors{iPred}).Corr = cat(1,GLM_coeff.summary.(Predictors{iPred}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Corr);
        GLM_coeff.summary.(Predictors{iPred}).Pvalue = cat(1,GLM_coeff.summary.(Predictors{iPred}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).ALL_Pvalue);
    for iPred2 = 1:3
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr(6:11,6:11));
        GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue = cat(2,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues(6:11,6:11));
    end
        GLM_coeff.summary.(Epoch{iEpoch}).Corr = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Corr,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Corr);
    GLM_coeff.summary.(Epoch{iEpoch}).Pvalue = cat(1,GLM_coeff.summary.(Epoch{iEpoch}).Pvalue,GLM_coeff.(Epoch{iEpoch}).(Predictors{iPred}).AllPred_Pvalue);
    end
end
%%
rColorMap = [linspace(233/255, 255/255, 183),linspace(255/255, 161/255, 73)];
    gColorMap = [linspace(163/255, 255/255, 183),linspace(255/255, 215/255, 73)];
    bColorMap = [linspace(201/255, 255/255, 183),linspace(255/255, 106/255, 73)];
colorMap = [rColorMap; gColorMap; bColorMap]';

mincolor = .01;
maxcolor = .2;
figure
for iPred = 1:3
    subplot(2,3,iPred)
    GLM_coeff.summary.(Predictors{iPred}).Pvalue(GLM_coeff.summary.(Predictors{iPred}).Pvalue == 0)=NaN;
    
heatmap(GLM_coeff.summary.(Predictors{iPred}).Pvalue,[],[],'%0.2','ColorMap', flipud(colorMap),...% 'Colorbar',true, ...
    'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0])%...
%    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
% set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
hold on
plot([6.5 6.5],[0.5 18.5],'k'); plot([.5 18.5],[6.5 6.5],'k');
plot([12.5 12.5],[0.5 18.5],'k'); plot([.5 18.5],[12.5 12.5],'k');
title(Predictors{iPred})
end

for iEpoch = 1:3
    subplot(2,3,iEpoch+3)
    GLM_coeff.summary.(Epoch{iEpoch}).Pvalue(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue == 0)=NaN;
    
heatmap(GLM_coeff.summary.(Epoch{iEpoch}).Pvalue,[],[],'%0.2','ColorMap', flipud(colorMap),...% 'Colorbar',true, ...
    'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0 0 0])%...
%    'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'});
% set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
hold on
plot([6.5 6.5],[0.5 18.5],'k'); plot([.5 18.5],[6.5 6.5],'k');
plot([12.5 12.5],[0.5 18.5],'k'); plot([.5 18.5],[12.5 12.5],'k');
title(Epoch{iEpoch})
end