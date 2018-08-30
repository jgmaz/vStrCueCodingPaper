function Class_accuracy_SHUFF = genLDAshuff(directory,destination,num_Shuffs,num_Iterations)
% function Class_accuracy_SHUFF = genLDAshuff(directory,destination,num_Shuffs,num_Iterations)
%
%
% INPUTS:
%
% OUTPUTS:

cd(directory)

%% get data into matrix for LDA
rng('shuffle')

mat_files = dir('*.mat');
t_count = 0;

for iTime = -.5:.1:.5
    t_count = t_count + 1;
    count = 1;
    cell_count = 1;
    
    for iShuff = 1:num_Shuffs
        LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(1:240,1:133) = NaN;
    end
    
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    for iCell = 1:length(mat_files)
        
        load(mat_files(iCell).name);
        
        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
            
            block_drift.block1_length(iCell) = length(FRATE.Cue.Trial_firing_rate_block1);
            block_drift.block1_half(iCell) = round(block_drift.block1_length(iCell) / 2);
            block_drift.b1_1st_avg(iCell) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(iCell)));
            block_drift.b1_2nd_avg(iCell) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(iCell)+1:end));
            block_drift.MWU_b1(iCell) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(iCell)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(iCell)+1:end));
            
            block_drift.block2_length(iCell) = length(FRATE.Cue.Trial_firing_rate_block2);
            block_drift.block2_half(iCell) = round(block_drift.block2_length(iCell) / 2);
            block_drift.b2_1st_avg(iCell) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(iCell)));
            block_drift.b2_2nd_avg(iCell) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(iCell)+1:end));
            block_drift.MWU_b2(iCell) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(iCell)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(iCell)+1:end));
            
            switch block_drift.MWU_b1(iCell) < .01 || block_drift.MWU_b2(iCell) < .01
                case 0
                    mat_overview.fname{iCell} = mat_files(iCell).name;
                    disp(cat(2,num2str(iCell),'/',num2str(length(dir('*.mat'))),' (',num2str(iTime),' epoch)'));
                    
                    
                    switch sesh.block_order %modality
                        case 1
                            starting_b1 = 0;
                            starting_b2 = 120;
                        case 2
                            starting_b1 = 120;
                            starting_b2 = 0;
                    end
                    
                    trials_count = 1;
                    
                    b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                    for jj = 1:b1_length
                        if metadata.TrialInfo{1,1}.trial_length_analysis(jj) < time_window_end
                            trial_length(trials_count) = metadata.TrialInfo{1,1}.trial_length_analysis(jj);
                        else
                            trial_length(trials_count) = time_window_end;
                        end
                        end_time(trials_count) = metadata.dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                        start_time(trials_count) = metadata.dataPoint.Trials(trials_count)/1000 + time_window_start;
                        % now, count spikes between start and end
                        these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                        % convert to firing rate and store
                        firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                        
                        trials_count = trials_count + 1;
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
                        end_time(trials_count) = metadata.dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
                        start_time(trials_count) = metadata.dataPoint.Trials(trials_count)/1000 + time_window_start;
                        % now, count spikes between start and end
                        these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                        % convert to firing rate and store
                        firing_rate(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
                        
                        trials_count = trials_count + 1;
                    end
                    for iShuff = 1:num_Shuffs
                        FRate_shuff = datasample(firing_rate,length(firing_rate),'Replace',false);
                        arm1_rew_count = 1;
                        arm2_rew_count = 16;
                        arm3_rew_count = 31;
                        arm4_rew_count = 46;
                        arm1_unrew_count = 61;
                        arm2_unrew_count = 76;
                        arm3_unrew_count = 91;
                        arm4_unrew_count = 106;
                        trials_count = 1;
                        
                        for jj = 1:b1_length
                            if jj > 100
                                break
                            end
                            
                            switch metadata.TrialInfo{1,1}.rewarded(jj)
                                case 0
                                    switch metadata.TrialInfo{1,1}.photosensorID(jj)
                                        case 1
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm1_unrew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm1_unrew_count = arm1_unrew_count + 1;
                                        case 2
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm2_unrew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm2_unrew_count = arm2_unrew_count + 1;
                                        case 3
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm3_unrew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm3_unrew_count = arm3_unrew_count + 1;
                                        case 4
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm4_unrew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm4_unrew_count = arm4_unrew_count + 1;
                                    end
                                case 1
                                    switch metadata.TrialInfo{1,1}.photosensorID(jj)
                                        case 1
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm1_rew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm1_rew_count = arm1_rew_count + 1;
                                        case 2
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm2_rew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm2_rew_count = arm2_rew_count + 1;
                                        case 3
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm3_rew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm3_rew_count = arm3_rew_count + 1;
                                        case 4
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm4_rew_count+starting_b1,cell_count) = FRate_shuff(trials_count);
                                            arm4_rew_count = arm4_rew_count + 1;
                                    end
                            end
                            
                            trials_count = trials_count + 1;
                            
                        end
                        
                        arm1_rew_count = 1;
                        arm2_rew_count = 16;
                        arm3_rew_count = 31;
                        arm4_rew_count = 46;
                        arm1_unrew_count = 61;
                        arm2_unrew_count = 76;
                        arm3_unrew_count = 91;
                        arm4_unrew_count = 106;
                        
                        for jj = 1:b2_length
                            if jj > 100
                                break
                            end
                            switch metadata.TrialInfo{1,2}.rewarded(jj)
                                case 0
                                    switch metadata.TrialInfo{1,2}.photosensorID(jj)
                                        case 1
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm1_unrew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm1_unrew_count = arm1_unrew_count + 1;
                                        case 2
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm2_unrew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm2_unrew_count = arm2_unrew_count + 1;
                                        case 3
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm3_unrew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm3_unrew_count = arm3_unrew_count + 1;
                                        case 4
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm4_unrew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm4_unrew_count = arm4_unrew_count + 1;
                                    end
                                case 1
                                    switch metadata.TrialInfo{1,2}.photosensorID(jj)
                                        case 1
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm1_rew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm1_rew_count = arm1_rew_count + 1;
                                        case 2
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm2_rew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm2_rew_count = arm2_rew_count + 1;
                                        case 3
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm3_rew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm3_rew_count = arm3_rew_count + 1;
                                        case 4
                                            LDA_raw_input{t_count}.(cat(2,'s',num2str(iShuff)))(arm4_rew_count+starting_b2,cell_count) = FRate_shuff(trials_count);
                                            arm4_rew_count = arm4_rew_count + 1;
                                    end
                            end
                            
                            trials_count = trials_count + 1;
                            
                        end
                    end
                    cell_count = cell_count + 1;
            end
        end
    end
end


%% 
mdl_identifier = {'Modality','Location','Outcome'};

for iLabel = 1:15
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 16:30
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 31:45
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 46:60
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 61:75
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 76:90
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 91:105
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 106:120
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 121:135
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 136:150
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 151:165
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 166:180
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
end
for iLabel = 181:195
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 196:210
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 211:225
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
end
for iLabel = 226:240
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
end

%%
for iTime = 1:length(LDA_raw_input)
    for iShuff = 1:num_Shuffs
        for iMdl = 1:length(mdl_identifier)
            disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
            Label_loop = unique(Labels.(mdl_identifier{iMdl}));
            Label_exclude = [];
            for iLabel = 1:length(Label_loop)
                Label_trials = strcmp(Labels.(mdl_identifier{iMdl}),Label_loop{iLabel});
                Label_variance = nanvar(LDA_raw_input{iTime}.(cat(2,'s',num2str(iShuff)))(Label_trials,:));
                Label_exclude = cat(2,Label_exclude,find(Label_variance < .15));
            end
            LDA_fixed_input.(mdl_identifier{iMdl}){iTime}.(cat(2,'s',num2str(iShuff))) = LDA_raw_input{iTime}.(cat(2,'s',num2str(iShuff)));
            LDA_fixed_input.(mdl_identifier{iMdl}){iTime}.(cat(2,'s',num2str(iShuff)))(:,Label_exclude) = [];
        end
    end
end

%%
for iTime = 1:length(LDA_raw_input) %note, skipped 4 and 5
    for iShuff = 1:num_Shuffs
        for iMdl = 1:length(mdl_identifier)
            disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
            for iteration = 1:num_Iterations
                LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}.(cat(2,'s',num2str(iShuff)))(1:240,:),Labels.(mdl_identifier{iMdl}));
                LDA_cv = crossval(LDA_Mdl);
                LDA_Error.(mdl_identifier{iMdl}).(cat(2,'s',num2str(iShuff)))(iteration,iTime) = kfoldLoss(LDA_cv);
            end
        end
    end
end

%%
for iShuff = 1:num_Shuffs
    for iMdl = 1:length(mdl_identifier)
        LDA_Error_SHUFF.(mdl_identifier{iMdl})(iShuff,:) = nanmean(LDA_Error.(mdl_identifier{iMdl}).(cat(2,'s',num2str(iShuff))));
    end
end

%%
for iMdl = 1:length(mdl_identifier)
    Class_accuracy_SHUFF.(mdl_identifier{iMdl}) = 1 - LDA_Error_SHUFF.(mdl_identifier{iMdl});
end

save(cat(2,destination,'LDA_SHUFF_performance.mat'),'Class_accuracy_SHUFF')

end