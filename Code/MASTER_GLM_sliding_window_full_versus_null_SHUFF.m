warning('error', 'stats:glmfit:IterationLimit'); % turn GLM warning into error so that you can use try catch to skip to next iteration of loop.
rng('shuffle')
%% time window to analyze
for iTime = -.5:.1:.5
    mdl_count = 1;
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    mat_files = dir('*.mat');
    for kk = 1:length(dir('*.mat'))
        load(mat_files(kk).name);
        mat_overview.fname{kk} = mat_files(kk).name;
        disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
        new_v_old = strcmp(mat_overview.fname{kk}(1:4),'R060');
        
        if RANK.two.Trial > 975 || RANK.two.Trial < 26
            if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                
                clear dataset ds %mdl
                switch new_v_old
                    case 0
                        %% old rats (R053,R056,R057)
                        trials_count = 1;
                        b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                        for jj = 1:b1_length
                            if metadata.TrialInfo_block1.trial_length_analysis(jj) < time_window_end
                                trial_length(trials_count) = metadata.TrialInfo_block1.trial_length_analysis(jj);
                            else
                                trial_length(trials_count) = time_window_end;
                            end
                            end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                            start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
                            
                            % now, count spikes between start and end
                            these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                            
                            % convert to firing rate and store
                            firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                            
                            %                         if metadata.TrialInfo_block1.summary(jj,2) == metadata.TrialInfo_block1.summary(jj,3) % if correct response
                            if jj == 1
                                dataset(trials_count,1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo_block1.rewarded(1)
                                    case 0
                                        dataset(trials_count,1) = 1;
                                    case 1
                                        dataset(trials_count,1) = 2;
                                end
                            else
                                switch metadata.TrialInfo_block1.rewarded(jj-1) + metadata.TrialInfo_block1.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count,1) = 0;
                                    case 1
                                        switch metadata.TrialInfo_block1.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count,1) = 1;
                                            case 1
                                                dataset(trials_count,1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count,1) = 3;
                                end
                            end
                            
                            dataset(trials_count,2) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                            switch sesh.block_order %modality
                                case 1
                                    dataset(trials_count,3) = 1;
                                case 2
                                    dataset(trials_count,3) = 2;
                            end
                            dataset(trials_count,4) = metadata.TrialInfo_block1.summary(jj,5); %arm
                            dataset(trials_count,5) = metadata.TrialInfo_block1.approached(jj); %behav - app
                            dataset(trials_count,6) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                            dataset(trials_count,7) = jj;
                            dataset(trials_count,8) = firing_rate(trials_count); %FR
                            
                            trials_count = trials_count + 1;
                            %                         end
                        end
                        
                        b2_start = b1_length + 1;
                        b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                        total = b1_length + b2_length;
                        
                        for jj = 1:b2_length
                            if metadata.TrialInfo_block2.trial_length_analysis(jj) < time_window_end
                                trial_length(trials_count) = metadata.TrialInfo_block2.trial_length_analysis(jj);
                            else
                                trial_length(trials_count) = time_window_end;
                            end
                            end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                            start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
                            
                            % now, count spikes between start and end
                            these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                            
                            % convert to firing rate and store
                            firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                            
                            %                         if metadata.TrialInfo_block2.summary(jj,2) == metadata.TrialInfo_block2.summary(jj,3)
                            if jj == 1
                                dataset(trials_count,1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo_block2.rewarded(1)
                                    case 0
                                        dataset(trials_count,1) = 1;
                                    case 1
                                        dataset(trials_count,1) = 2;
                                end
                            else
                                switch metadata.TrialInfo_block2.rewarded(jj-1) + metadata.TrialInfo_block2.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count,1) = 0;
                                    case 1
                                        switch metadata.TrialInfo_block2.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count,1) = 1;
                                            case 1
                                                dataset(trials_count,1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count,1) = 3;
                                end
                            end
                            
                            dataset(trials_count,2) = metadata.TrialInfo_block2.rewarded(jj);
                            switch sesh.block_order
                                case 1
                                    dataset(trials_count,3) = 2;
                                case 2
                                    dataset(trials_count,3) = 1;
                            end
                            dataset(trials_count,4) = metadata.TrialInfo_block2.summary(jj,5);
                            dataset(trials_count,5) = metadata.TrialInfo_block2.approached(jj);
                            dataset(trials_count,6) = metadata.TrialInfo_block2.summary(jj,15);
                            dataset(trials_count,7) = b1_length+jj;
                            dataset(trials_count,8) = firing_rate(trials_count);
                            
                            trials_count = trials_count + 1;
                            %                         end
                        end
                        
                    case 1
                        %% new rats (R060)
                        trials_count = 1;
                        b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                        for jj = 1:b1_length
                            if metadata.TrialInfo{1,1}.trial_length_analysis(jj) < time_window_end
                                trial_length(trials_count) = metadata.TrialInfo{1,1}.trial_length_analysis(jj);
                            else
                                trial_length(trials_count) = time_window_end;
                            end
                            end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                            start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
                            
                            % now, count spikes between start and end
                            these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                            
                            % convert to firing rate and store
                            firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                            
                            %                         if metadata.TrialInfo{1,1}.summary(jj,2) == metadata.TrialInfo{1,1}.summary(jj,3)
                            if jj == 1
                                dataset(trials_count,1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo{1,1}.rewarded(1)
                                    case 0
                                        dataset(trials_count,1) = 1;
                                    case 1
                                        dataset(trials_count,1) = 2;
                                end
                            else
                                switch metadata.TrialInfo{1,1}.rewarded(jj-1) + metadata.TrialInfo{1,1}.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count,1) = 0;
                                    case 1
                                        switch metadata.TrialInfo{1,1}.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count,1) = 1;
                                            case 1
                                                dataset(trials_count,1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count,1) = 3;
                                end
                            end
                            
                            dataset(trials_count,2) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                            switch sesh.block_order %modality
                                case 1
                                    dataset(trials_count,3) = 1;
                                case 2
                                    dataset(trials_count,3) = 2;
                            end
                            dataset(trials_count,4) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                            dataset(trials_count,5) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                            dataset(trials_count,6) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                            dataset(trials_count,7) = jj;
                            dataset(trials_count,8) = firing_rate(trials_count); %FR
                            
                            trials_count = trials_count + 1;
                            %                         end
                        end
                        
                        b2_start = b1_length + 1;
                        b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                        total = b1_length + b2_length;
                        
                        for jj = 1:b2_length
                            if metadata.TrialInfo{1,2}.trial_length_analysis(jj) < time_window_end
                                trial_length(trials_count) = metadata.TrialInfo{1,2}.trial_length_analysis(jj);
                            else
                                trial_length(trials_count) = time_window_end;
                            end
                            end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                            start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
                            
                            % now, count spikes between start and end
                            these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                            
                            % convert to firing rate and store
                            firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                            
                            %                         if metadata.TrialInfo{1,2}.summary(jj,2) == metadata.TrialInfo{1,2}.summary(jj,3)
                            if jj == 1
                                dataset(trials_count,1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo{1,2}.rewarded(1)
                                    case 0
                                        dataset(trials_count,1) = 1;
                                    case 1
                                        dataset(trials_count,1) = 2;
                                end
                            else
                                switch metadata.TrialInfo{1,2}.rewarded(jj-1) + metadata.TrialInfo{1,2}.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count,1) = 0;
                                    case 1
                                        switch metadata.TrialInfo{1,2}.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count,1) = 1;
                                            case 1
                                                dataset(trials_count,1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count,1) = 3;
                                end
                            end
                            
                            dataset(trials_count,2) = metadata.TrialInfo{1,2}.rewarded(jj);
                            switch sesh.block_order
                                case 1
                                    dataset(trials_count,3) = 2;
                                case 2
                                    dataset(trials_count,3) = 1;
                            end
                            dataset(trials_count,4) = metadata.TrialInfo{1,2}.summary(jj,5);
                            dataset(trials_count,5) = metadata.TrialInfo{1,2}.approached(jj);
                            dataset(trials_count,6) = metadata.TrialInfo{1,2}.summary(jj,15);
                            dataset(trials_count,7) = b1_length+jj;
                            dataset(trials_count,8) = firing_rate(trials_count);
                            
                            trials_count = trials_count + 1;
                            %                         end
                        end
                        
                end
                
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
        
        if block_drift.MWU_b1(kk) > .01 && block_drift.MWU_b2(kk) > .01
                
                if sum(dataset(:,8)) > 10

                        for iShuff = 1:25
                            dataset2 = dataset;
                           dataset2(:,8) = datasample(dataset(:,8),length(dataset(:,8)),'Replace',false);                    
                    
                    ds = array2table(dataset2,'VariableNames',{'Previous','Outcome','Modality','Location','Approach','Latency','Trial','FiringRate'});
                    modelspec = {'FiringRate ~ 1 + Outcome + Modality + Location + Previous + Approach + Latency + Trial',...
                        'FiringRate ~ 1 + Outcome + Location + Previous + Approach + Latency + Trial',...
                        'FiringRate ~ 1 + Outcome + Modality + Previous + Approach + Latency + Trial',...
                        'FiringRate ~ 1 + Modality + Location + Previous + Approach + Latency + Trial'};
                    mdlidentifier = {'Full','Modality','Location','Outcome'};
                    mdl_cell{mdl_count} = mat_overview.fname{kk};
                    for iMdl = 1:length(mdlidentifier)
                        %             mdl{kk}= stepwiseglm(ds,'constant','upper','linear','Distribution','poisson');
                        try
                            mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count} = fitglm(ds,modelspec{iMdl},'Distribution','poisson');
                            mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl})(mdl_count,1) = (mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{1}){mdl_count}.Rsquared.Adjusted - mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count}.Rsquared.Adjusted) * 100;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl})(mdl_count,1) = mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{1}){mdl_count}.ModelCriterion.AIC - mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count}.ModelCriterion.AIC;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl})(mdl_count,1) = mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{1}){mdl_count}.ModelCriterion.BIC - mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count}.ModelCriterion.BIC;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl})(mdl_count,1) = mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{1}){mdl_count}.Deviance - mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count}.Deviance;
                        catch %err
                            mdl.(cat(2,'s',num2str(iShuff))).(mdlidentifier{iMdl}){mdl_count} = [];
                            mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl})(mdl_count,1) = NaN;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl})(mdl_count,1) = NaN;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl})(mdl_count,1) = NaN;
                            mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl})(mdl_count,1) = NaN;
                        end
                    end
                        end
                    mdl_count = mdl_count + 1;
                end
                
        end
            end
        end
    end
save(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-06-28-GLM_cueon_sliding_window_full_v_null_',num2str(iTime),'_SHUFF.mat'),'mdl','mdlcomp')
clearvars -except iTime
end

%%
t_count = 0;
mdlidentifier = {'Full','Modality','Location','Outcome'};
for iTime = -.5:.1:.5
    load(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-06-28-GLM_cueon_sliding_window_full_v_null_',num2str(iTime),'_SHUFF.mat'));
    t_count = t_count + 1;
    disp(num2str(iTime))
    for iShuff = 1:25
    for iMdl = 1:3
    mdlavg.(cat(2,'s',num2str(iShuff))).Rsquared.MEAN(t_count,iMdl) = nanmean(mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl+1}));
    mdlavg.(cat(2,'s',num2str(iShuff))).Rsquared.SEM(t_count,iMdl) = nanstd(mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.(cat(2,'s',num2str(iShuff))).Rsquared.(mdlidentifier{iMdl+1}))));
    mdlavg.(cat(2,'s',num2str(iShuff))).AIC.MEAN(t_count,iMdl) = nanmean(mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl+1}));
    mdlavg.(cat(2,'s',num2str(iShuff))).AIC.SEM(t_count,iMdl) = nanstd(mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.(cat(2,'s',num2str(iShuff))).AIC.(mdlidentifier{iMdl+1}))));
    mdlavg.(cat(2,'s',num2str(iShuff))).BIC.MEAN(t_count,iMdl) = nanmean(mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl+1}));
    mdlavg.(cat(2,'s',num2str(iShuff))).BIC.SEM(t_count,iMdl) = nanstd(mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.(cat(2,'s',num2str(iShuff))).BIC.(mdlidentifier{iMdl+1}))));
    mdlavg.(cat(2,'s',num2str(iShuff))).Deviance.MEAN(t_count,iMdl) = nanmean(mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl+1}));
    mdlavg.(cat(2,'s',num2str(iShuff))).Deviance.SEM(t_count,iMdl) = nanstd(mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.(cat(2,'s',num2str(iShuff))).Deviance.(mdlidentifier{iMdl+1}))));
    end
    end
end

%%
t_count = 0;
mdlidentifier = {'Modality','Location','Outcome'};
mdlcriterion = {'Rsquared','AIC','BIC','Deviance'};
for iTime = -.5:.1:.5
        t_count = t_count + 1;
        for iShuff = 1:25
            for iMdl = 1:3
                for iCrit = 1:4
        mdlavgSHUFF.(mdlcriterion{iCrit}).(mdlidentifier{iMdl})(t_count,iShuff) = mdlavg.(cat(2,'s',num2str(iShuff))).(mdlcriterion{iCrit}).MEAN(t_count,iMdl);
                end
            end
        end
end

%%
for iTime = 1:length(mdlavgSHUFF.(mdlcriterion{1}).(mdlidentifier{1})(:,1))
for iMdl = 1:3
    for iCrit = 1:4
        mdlavgSHUFF.(mdlcriterion{iCrit}).MEAN(iTime,iMdl) = mean(mdlavgSHUFF.(mdlcriterion{iCrit}).(mdlidentifier{iMdl})(iTime,:));
        mdlavgSHUFF.(mdlcriterion{iCrit}).SEM(iTime,iMdl) = std(mdlavgSHUFF.(mdlcriterion{iCrit}).(mdlidentifier{iMdl})(iTime,:))/sqrt(length(mdlavgSHUFF.(mdlcriterion{1}).(mdlidentifier{1})));
    end
end
end

%%
t_count = 0;
mdlidentifier = {'Full','Modality','Location','Outcome'};
for iTime = -.5:.1:.5
    load(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-05-25-GLM_cueon_sliding_window_full_v_null_',num2str(iTime),'.mat'));
    t_count = t_count + 1;
    for iMdl = 1:3
    mdlavg.Rsquared.MEAN(t_count,iMdl) = nanmean(mdlcomp.Rsquared.(mdlidentifier{iMdl+1}));
    mdlavg.Rsquared.SEM(t_count,iMdl) = nanstd(mdlcomp.Rsquared.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.Rsquared.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.Rsquared.(mdlidentifier{iMdl+1}))));
    mdlavg.AIC.MEAN(t_count,iMdl) = nanmean(mdlcomp.AIC.(mdlidentifier{iMdl+1}));
    mdlavg.AIC.SEM(t_count,iMdl) = nanstd(mdlcomp.AIC.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.AIC.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.AIC.(mdlidentifier{iMdl+1}))));
    mdlavg.BIC.MEAN(t_count,iMdl) = nanmean(mdlcomp.BIC.(mdlidentifier{iMdl+1}));
    mdlavg.BIC.SEM(t_count,iMdl) = nanstd(mdlcomp.BIC.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.BIC.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.BIC.(mdlidentifier{iMdl+1}))));
    mdlavg.Deviance.MEAN(t_count,iMdl) = nanmean(mdlcomp.Deviance.(mdlidentifier{iMdl+1}));
    mdlavg.Deviance.SEM(t_count,iMdl) = nanstd(mdlcomp.Deviance.(mdlidentifier{iMdl+1})) / sqrt(numel(mdlcomp.Deviance.(mdlidentifier{iMdl+1}))-sum(isnan(mdlcomp.Deviance.(mdlidentifier{iMdl+1}))));
    end
end

%% 
figure;
subplot(2,2,1)
shadedErrorBar(-.25:.1:.75,mdlavg.Rsquared.MEAN(:,1) ,mdlavg.Rsquared.SEM(:,1),'-b',1);
hold on
shadedErrorBar(-.25:.1:.75,mdlavg.Rsquared.MEAN(:,2) ,mdlavg.Rsquared.SEM(:,2),'-r',1);
shadedErrorBar(-.25:.1:.75,mdlavg.Rsquared.MEAN(:,3) ,mdlavg.Rsquared.SEM(:,3),'-y',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Rsquared.MEAN(:,1) ,mdlavgSHUFF.Rsquared.SEM(:,1),'--b',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Rsquared.MEAN(:,2) ,mdlavgSHUFF.Rsquared.SEM(:,2),'--r',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Rsquared.MEAN(:,3) ,mdlavgSHUFF.Rsquared.SEM(:,3),'--y',1);
plot(.2,0.1:.1:4,'.k'); plot(-.2,0.1:.1:4,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('Rsquared'); ylabel('Average improvement to Rsquared')
% ylim([0 .5]);
xlabel('Time from cue-onset (s)');

subplot(2,2,2)
shadedErrorBar(-.25:.1:.75,mdlavg.AIC.MEAN(:,1) ,mdlavg.AIC.SEM(:,1),'-b',1);
hold on
shadedErrorBar(-.25:.1:.75,mdlavg.AIC.MEAN(:,2) ,mdlavg.AIC.SEM(:,2),'-r',1);
shadedErrorBar(-.25:.1:.75,mdlavg.AIC.MEAN(:,3) ,mdlavg.AIC.SEM(:,3),'-y',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.AIC.MEAN(:,1) ,mdlavgSHUFF.AIC.SEM(:,1),'--b',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.AIC.MEAN(:,2) ,mdlavgSHUFF.AIC.SEM(:,2),'--r',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.AIC.MEAN(:,3) ,mdlavgSHUFF.AIC.SEM(:,3),'--y',1);
plot(.2,0:-1:-45,'.k'); plot(-.2,0:-1:-45,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('AIC'); ylabel('Average improvement to AIC')
xlabel('Time from cue-onset (s)');

subplot(2,2,3)
shadedErrorBar(-.25:.1:.75,mdlavg.BIC.MEAN(:,1) ,mdlavg.BIC.SEM(:,1),'-b',1);
hold on
shadedErrorBar(-.25:.1:.75,mdlavg.BIC.MEAN(:,2) ,mdlavg.BIC.SEM(:,2),'-r',1);
shadedErrorBar(-.25:.1:.75,mdlavg.BIC.MEAN(:,3) ,mdlavg.BIC.SEM(:,3),'-y',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.BIC.MEAN(:,1) ,mdlavgSHUFF.BIC.SEM(:,1),'--b',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.BIC.MEAN(:,2) ,mdlavgSHUFF.BIC.SEM(:,2),'--r',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.BIC.MEAN(:,3) ,mdlavgSHUFF.BIC.SEM(:,3),'--y',1);
plot(.2,0:-1:-40,'.k'); plot(-.2,0:-1:-40,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('BIC'); ylabel('Average improvement to BIC')
xlabel('Time from cue-onset (s)');

subplot(2,2,4)
shadedErrorBar(-.25:.1:.75,mdlavg.Deviance.MEAN(:,1) ,mdlavg.Deviance.SEM(:,1),'-b',1);
hold on
shadedErrorBar(-.25:.1:.75,mdlavg.Deviance.MEAN(:,2) ,mdlavg.Deviance.SEM(:,2),'-r',1);
shadedErrorBar(-.25:.1:.75,mdlavg.Deviance.MEAN(:,3) ,mdlavg.Deviance.SEM(:,3),'-y',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Deviance.MEAN(:,1) ,mdlavgSHUFF.Deviance.SEM(:,1),'--b',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Deviance.MEAN(:,2) ,mdlavgSHUFF.Deviance.SEM(:,2),'--r',1);
shadedErrorBar(-.25:.1:.75,mdlavgSHUFF.Deviance.MEAN(:,3) ,mdlavgSHUFF.Deviance.SEM(:,3),'--y',1);
plot(.2,0:-1:-45,'.k'); plot(-.2,0:-1:-45,'.k');
%  legend({'identity' 'location' 'outcome'}); 
title('Deviance'); ylabel('Average improvement to Deviance')
xlabel('Time from cue-onset (s)');