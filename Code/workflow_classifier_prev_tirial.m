%% get data into matrix for LDA


mat_files = dir('*.mat');
t_count = 0;

for iTime = -.5:.1:.5
    t_count = t_count + 1;
    count = 1;
    cell_count = 1;
    
    LDA_raw_input{t_count}(1:240,1:133) = NaN;
    
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    for iCell = 1:length(mat_files)
        
        load(mat_files(iCell).name);
        
        if RANK.two.Trial > 975 || RANK.two.Trial < 26
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
                        new_v_old = strcmp(mat_overview.fname{iCell}(1:4),'R060');
                        
                        
                        switch sesh.block_order %modality
                            case 1
                                starting_b1 = 0;
                                starting_b2 = 120;
                            case 2
                                starting_b1 = 120;
                                starting_b2 = 0;
                        end
                        
                        switch new_v_old
                            case 0
                                %% old rats (R053,R056,R057)
                                trials_count = 1;
                                arm1_rew_count = 1;
                                arm2_rew_count = 16;
                                arm3_rew_count = 31;
                                arm4_rew_count = 46;
                                arm1_unrew_count = 61;
                                arm2_unrew_count = 76;
                                arm3_unrew_count = 91;
                                arm4_unrew_count = 106;
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    if jj > 100
                                        break
                                    end
%                                         if metadata.TrialInfo_block1.trial_length_analysis(jj) < time_window_end
%                                             trial_length(trials_count) = metadata.TrialInfo_block1.trial_length_analysis(jj);
%                                         else
%                                             trial_length(trials_count) = time_window_end;
%                                         end
%                                         end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
%                                         start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
%                                         % now, count spikes between start and end
%                                         these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
if jj == 1
    prev_trial(trials_count) = NaN;
else
                                        % convert to firing rate and store
                                        prev_trial(trials_count) = metadata.TrialInfo_block1.rewarded(jj-1);
end
                                        
                                        switch metadata.TrialInfo_block1.rewarded(jj)
                                            case 0
                                                switch metadata.TrialInfo_block1.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm1_unrew_count = arm1_unrew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm2_unrew_count = arm2_unrew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm3_unrew_count = arm3_unrew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm4_unrew_count = arm4_unrew_count + 1;
                                                end
                                            case 1
                                                switch metadata.TrialInfo_block1.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm1_rew_count = arm1_rew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm2_rew_count = arm2_rew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm3_rew_count = arm3_rew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm4_rew_count = arm4_rew_count + 1;
                                                end
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
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
%                                         if metadata.TrialInfo_block2.trial_length_analysis(jj) < time_window_end
%                                             trial_length(trials_count) = metadata.TrialInfo_block2.trial_length_analysis(jj);
%                                         else
%                                             trial_length(trials_count) = time_window_end;
%                                         end
%                                         end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
%                                         start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
%                                         % now, count spikes between start and end
%                                         these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                                        % convert to firing rate and store
                                        if jj == 1
                                            prev_trial(trials_count) = NaN;
                                        else
                                        prev_trial(trials_count) = metadata.TrialInfo_block2.rewarded(jj - 1);
                                        end
                                        
                                        switch metadata.TrialInfo_block2.rewarded(jj)
                                            case 0
                                                switch metadata.TrialInfo_block2.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm1_unrew_count = arm1_unrew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm2_unrew_count = arm2_unrew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm3_unrew_count = arm3_unrew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm4_unrew_count = arm4_unrew_count + 1;
                                                end
                                            case 1
                                                switch metadata.TrialInfo_block2.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm1_rew_count = arm1_rew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm2_rew_count = arm2_rew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm3_rew_count = arm3_rew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm4_rew_count = arm4_rew_count + 1;
                                                end
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                            case 1
                                %% new rats (R060)
                                trials_count = 1;
                                arm1_rew_count = 1;
                                arm2_rew_count = 16;
                                arm3_rew_count = 31;
                                arm4_rew_count = 46;
                                arm1_unrew_count = 61;
                                arm2_unrew_count = 76;
                                arm3_unrew_count = 91;
                                arm4_unrew_count = 106;
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    if jj > 100
                                        break
                                    end
%                                     if metadata.TrialInfo{1,1}.trial_length_analysis(jj) < time_window_end
%                                         trial_length(trials_count) = metadata.TrialInfo{1,1}.trial_length_analysis(jj);
%                                     else
%                                         trial_length(trials_count) = time_window_end;
%                                     end
%                                     end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
%                                     start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
%                                     % now, count spikes between start and end
%                                     these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
%                                     % convert to firing rate and store
%                                     prev_trial(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
%                                     
if jj == 1
                                            prev_trial(trials_count) = NaN;
else
                                        prev_trial(trials_count) = metadata.TrialInfo{1,1}.rewarded(jj - 1);
end
                                        
                                    switch metadata.TrialInfo{1,1}.rewarded(jj)
                                        case 0
                                            switch metadata.TrialInfo{1,1}.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm1_unrew_count = arm1_unrew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm2_unrew_count = arm2_unrew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm3_unrew_count = arm3_unrew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_unrew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm4_unrew_count = arm4_unrew_count + 1;
                                                end
                                        case 1
                                            switch metadata.TrialInfo{1,1}.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm1_rew_count = arm1_rew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm2_rew_count = arm2_rew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm3_rew_count = arm3_rew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_rew_count+starting_b1,cell_count) = prev_trial(trials_count);
                                                        arm4_rew_count = arm4_rew_count + 1;
                                                end
                                    end
                                    
                                    trials_count = trials_count + 1;
                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
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
%                                     if metadata.TrialInfo{1,2}.trial_length_analysis(jj) < time_window_end
%                                         trial_length(trials_count) = metadata.TrialInfo{1,2}.trial_length_analysis(jj);
%                                     else
%                                         trial_length(trials_count) = time_window_end;
%                                     end
%                                     end_time(trials_count) = dataPoint.Trials(trials_count)/1000 + trial_length(trials_count);
%                                     start_time(trials_count) = dataPoint.Trials(trials_count)/1000 + time_window_start;
%                                     % now, count spikes between start and end
%                                     these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
%                                     % convert to firing rate and store
%                                     prev_trial(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
if jj == 1
                                            prev_trial(trials_count) = NaN;
else
                                        prev_trial(trials_count) = metadata.TrialInfo{1,2}.rewarded(jj - 1);                                    
end
                                        
                                    switch metadata.TrialInfo{1,2}.rewarded(jj)
                                        case 0
                                            switch metadata.TrialInfo{1,2}.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm1_unrew_count = arm1_unrew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm2_unrew_count = arm2_unrew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm3_unrew_count = arm3_unrew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_unrew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm4_unrew_count = arm4_unrew_count + 1;
                                                end
                                        case 1
                                            switch metadata.TrialInfo{1,2}.photosensorID(jj)
                                                    case 1
                                                        LDA_raw_input{t_count}(arm1_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm1_rew_count = arm1_rew_count + 1;
                                                    case 2
                                                        LDA_raw_input{t_count}(arm2_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm2_rew_count = arm2_rew_count + 1;
                                                    case 3
                                                        LDA_raw_input{t_count}(arm3_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
                                                        arm3_rew_count = arm3_rew_count + 1;
                                                    case 4
                                                        LDA_raw_input{t_count}(arm4_rew_count+starting_b2,cell_count) = prev_trial(trials_count);
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
end


%% Figure: mod; loc; out; modxloc; modxout; locxout; modxlocxout;
mdl_identifier = {'Modality','Location','Outcome','ModxOut','ModxLoc','LocxOut','ModxLocxOut'};

for iLabel = 1:15    
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'L1';
    Labels.LocxOut{iLabel,1} = '1+';
    Labels.ModxLocxOut{iLabel,1} = 'L1+';
    Labels.ModxOut{iLabel,1} = 'L+';
    if iLabel < 7
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 6 && iLabel < 13
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 16:30
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'L2';
    Labels.LocxOut{iLabel,1} = '2+';
    Labels.ModxLocxOut{iLabel,1} = 'L2+';
    Labels.ModxOut{iLabel,1} = 'L+';
    if iLabel < 22
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 21 && iLabel < 28
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 31:45
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'L3';
    Labels.LocxOut{iLabel,1} = '3+';
    Labels.ModxLocxOut{iLabel,1} = 'L3+';
    Labels.ModxOut{iLabel,1} = 'L+';
    if iLabel < 37
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 36 && iLabel < 43
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 46:60
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'L4';
    Labels.LocxOut{iLabel,1} = '4+';
    Labels.ModxLocxOut{iLabel,1} = 'L4+';
    Labels.ModxOut{iLabel,1} = 'L+';
    if iLabel < 52
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 51 && iLabel < 58
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 61:75
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'L1';
    Labels.LocxOut{iLabel,1} = '1-';
    Labels.ModxLocxOut{iLabel,1} = 'L1-';
    Labels.ModxOut{iLabel,1} = 'L-';
    if iLabel < 67
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 66 && iLabel < 73
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 76:90
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'L2';
    Labels.LocxOut{iLabel,1} = '2-';
    Labels.ModxLocxOut{iLabel,1} = 'L2-';
    Labels.ModxOut{iLabel,1} = 'L-';
    if iLabel < 82
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 81 && iLabel < 88
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 91:105
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'L3';
    Labels.LocxOut{iLabel,1} = '3-';
    Labels.ModxLocxOut{iLabel,1} = 'L3-';
    Labels.ModxOut{iLabel,1} = 'L-';
    if iLabel < 97
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 96 && iLabel < 103
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 106:120
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'L';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'L4';
    Labels.LocxOut{iLabel,1} = '4-';
    Labels.ModxLocxOut{iLabel,1} = 'L4-';
    Labels.ModxOut{iLabel,1} = 'L-';
    if iLabel < 112
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 111 && iLabel < 118
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 121:135
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'S1';
    Labels.LocxOut{iLabel,1} = '1+';
    Labels.ModxLocxOut{iLabel,1} = 'S1+';
    Labels.ModxOut{iLabel,1} = 'S+';
    if iLabel < 127
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 126 && iLabel < 133
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 136:150
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'S2';
    Labels.LocxOut{iLabel,1} = '2+';
    Labels.ModxLocxOut{iLabel,1} = 'S2+';
    Labels.ModxOut{iLabel,1} = 'S+';
    if iLabel < 142
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 141 && iLabel < 148
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 151:165
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'S3';
    Labels.LocxOut{iLabel,1} = '3+';
    Labels.ModxLocxOut{iLabel,1} = 'S3+';
    Labels.ModxOut{iLabel,1} = 'S+';
    if iLabel < 157
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 156 && iLabel < 163
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 166:180
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '+';
    Labels.ModxLoc{iLabel,1} = 'S4';
    Labels.LocxOut{iLabel,1} = '4+';
    Labels.ModxLocxOut{iLabel,1} = 'S4+';
    Labels.ModxOut{iLabel,1} = 'S+';
    if iLabel < 172
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 171 && iLabel < 178
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 181:195
    Labels.Location{iLabel,1} = 'Arm 1';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'S1';
    Labels.LocxOut{iLabel,1} = '1-';
    Labels.ModxLocxOut{iLabel,1} = 'S1-';
    Labels.ModxOut{iLabel,1} = 'S-';
    if iLabel < 187
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 186 && iLabel < 193
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 196:210
    Labels.Location{iLabel,1} = 'Arm 2';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'S2';
    Labels.LocxOut{iLabel,1} = '2-';
    Labels.ModxLocxOut{iLabel,1} = 'S2-';
    Labels.ModxOut{iLabel,1} = 'S-';
    if iLabel < 202
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 201 && iLabel < 208
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 211:225
    Labels.Location{iLabel,1} = 'Arm 3';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'S3';
    Labels.LocxOut{iLabel,1} = '3-';
    Labels.ModxLocxOut{iLabel,1} = 'S3-';
    Labels.ModxOut{iLabel,1} = 'S-';
    if iLabel < 217
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 216 && iLabel < 223
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end
for iLabel = 226:240
    Labels.Location{iLabel,1} = 'Arm 4';
    Labels.Modality{iLabel,1} = 'S';
    Labels.Outcome{iLabel,1} = '-';
    Labels.ModxLoc{iLabel,1} = 'S4';
    Labels.LocxOut{iLabel,1} = '4-';
    Labels.ModxLocxOut{iLabel,1} = 'S4-';
    Labels.ModxOut{iLabel,1} = 'S-';
    if iLabel < 232
        Labels.Drift{iLabel,1} = '1st';
    elseif iLabel > 231 && iLabel < 238
        Labels.Drift{iLabel,1} = '2nd';
    else
        Labels.Drift{iLabel,1} = '2nd';
    end
end

%% remove cells with 0 variance for each label
for iTime = 1:length(LDA_raw_input) %note, skipped 4 and 5
    for iMdl = 1:3%length(mdl_identifier)
        disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
        Label_loop = unique(Labels.(mdl_identifier{iMdl}));
        Label_exclude = [];
        for iLabel = 1:length(Label_loop)
           Label_trials = strcmp(Labels.(mdl_identifier{iMdl}),Label_loop{iLabel});
               Label_variance = nanvar(LDA_raw_input{iTime}(Label_trials,:));
               Label_exclude = cat(2,Label_exclude,find(Label_variance < .15));
        end
        LDA_fixed_input.(mdl_identifier{iMdl}){iTime} = LDA_raw_input{iTime};
        LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(:,Label_exclude) = [];       
    end
end
       
%% Cross-validation
rng('shuffle')

for iTime = 1:length(LDA_raw_input) %note, skipped 4 and 5
    for iMdl = 1:3%length(mdl_identifier)
        disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))        
    for iteration = 1:100
    LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(1:240,:),Labels.(mdl_identifier{iMdl}));
    LDA_cv = crossval(LDA_Mdl);
    LDA_Error.(mdl_identifier{iMdl})(iteration,iTime) = kfoldLoss(LDA_cv);
    end   
    end
end

% [1:91 93:96 98:124 126:133]
%% Convert MCR to CR
for iMdl = 1:3%length(mdl_identifier)
Class_accuracy.(mdl_identifier{iMdl}) = 1 - LDA_Error.(mdl_identifier{iMdl});
end

%% Plot LDA
figure
for iMdl = 1:3%length(mdl_identifier)
subplot(3,3,iMdl)
violin(Class_accuracy.(mdl_identifier{iMdl}));
set(gca,'XTickLabel',{'','','','','','','','','','',''})
ylim([0 1])
hold on;
switch iMdl
    case 1
plot(0:.25:22,.5,'.','color','black');
    case 2
plot(0:.25:22,.25,'.','color','black');
    case 3
plot(0:.25:22,.5,'.','color','black');
    case 4
plot(0:.25:22,.125,'.','color','black');
ylabel('Classification accuracy')
    case 5
plot(0:.25:22,.25,'.','color','black');
    case 6
plot(0:.25:22,.125,'.','color','black');
    case 7
plot(0:.25:22,.0625,'.','color','black');
xlabel('Starting time')
end
title(mdl_identifier{iMdl})
set(gca,'FontSize',18);
end