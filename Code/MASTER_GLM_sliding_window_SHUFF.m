warning('error', 'stats:glmfit:IterationLimit'); % turn GLM warning into error so that you can use try catch to skip to next iteration of loop.
rng('shuffle')
%% time window to analyze
    for iTime = -.4:.1:.5


time_window_start = iTime; %starting time window for analysis, 0 = time zero
        time_window_end = time_window_start + .5;
        epoch_start = 0;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
%     if kk == 98 || kk == 4;
%         mdl{kk} = [];
%         continue
%     end
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
                    trials_count.(cat(2,'shuff_',num2str(iShuff))) =  1;
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
                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo_block1.rewarded(1)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                    case 1
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                end
                            else
                                switch metadata.TrialInfo_block1.rewarded(jj-1) + metadata.TrialInfo_block1.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 0;
                                    case 1
                                        switch metadata.TrialInfo_block1.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                            case 1
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 3;
                                end
                            end
                            
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),2) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                            switch sesh.block_order %modality
                                case 1
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                case 2
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 2;
                            end
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),4) = metadata.TrialInfo_block1.summary(jj,5); %arm
                                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),5) = metadata.TrialInfo_block1.approached(jj); %behav - app
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),6) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),7) = jj;
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),8) = firing_rate(trials_count); %FR
                            
                            trials_count.(cat(2,'shuff_',num2str(iShuff))) =  trials_count.(cat(2,'shuff_',num2str(iShuff))) + 1;
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
                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo_block2.rewarded(1)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                    case 1
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                end
                            else
                                switch metadata.TrialInfo_block2.rewarded(jj-1) + metadata.TrialInfo_block2.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 0;
                                    case 1
                                        switch metadata.TrialInfo_block2.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                            case 1
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 3;
                                end
                            end
                            
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),2) = metadata.TrialInfo_block2.rewarded(jj);
                            switch sesh.block_order
                                case 1
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 2;
                                case 2
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                            end
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),4) = metadata.TrialInfo_block2.summary(jj,5);
                                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),5) = metadata.TrialInfo_block2.approached(jj);
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),6) = metadata.TrialInfo_block2.summary(jj,15);
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),7) = b1_length+jj;
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),8) = firing_rate(trials_count);
                            
                            trials_count.(cat(2,'shuff_',num2str(iShuff))) =  trials_count.(cat(2,'shuff_',num2str(iShuff))) + 1;
%                         end
                    end
                    
                case 1
                    %% new rats (R060)
                    trials_count.(cat(2,'shuff_',num2str(iShuff))) =  1;
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
                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo{1,1}.rewarded(1)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                    case 1
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                end
                            else
                                switch metadata.TrialInfo{1,1}.rewarded(jj-1) + metadata.TrialInfo{1,1}.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 0;
                                    case 1
                                        switch metadata.TrialInfo{1,1}.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                            case 1
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 3;
                                end
                            end
                            
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),2) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                            switch sesh.block_order %modality
                                case 1
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                case 2
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 2;
                            end
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),4) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                                                         dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),5) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),6) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),7) = jj;
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),8) = firing_rate(trials_count); %FR
                            
                            trials_count.(cat(2,'shuff_',num2str(iShuff))) =  trials_count.(cat(2,'shuff_',num2str(iShuff))) + 1;
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
                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = NaN; %prev trial
                            elseif jj == 2
                                switch metadata.TrialInfo{1,2}.rewarded(1)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                    case 1
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                end
                            else
                                switch metadata.TrialInfo{1,2}.rewarded(jj-1) + metadata.TrialInfo{1,2}.rewarded(jj-2)
                                    case 0
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 0;
                                    case 1
                                        switch metadata.TrialInfo{1,2}.rewarded(jj-1)
                                            case 0
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                            case 1
                                                dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 2;
                                        end
                                    case 2
                                        dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),1) = 3;
                                end
                            end
                            
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),2) = metadata.TrialInfo{1,2}.rewarded(jj);
                            switch sesh.block_order
                                case 1
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 2;
                                case 2
                                    dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                            end
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),4) = metadata.TrialInfo{1,2}.summary(jj,5);
                                                         dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),5) = metadata.TrialInfo{1,2}.approached(jj);
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),6) = metadata.TrialInfo{1,2}.summary(jj,15);
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),7) = b1_length+jj;
                            dataset(trials_count.(cat(2,'shuff_',num2str(iShuff))),8) = firing_rate(trials_count);
                            
                            trials_count.(cat(2,'shuff_',num2str(iShuff))) =  trials_count.(cat(2,'shuff_',num2str(iShuff))) + 1;
%                         end
                    end
                    
            end
            if sum(dataset(:,8)) < 30
                mdl{kk} = [];
            else
                %%
                for iShuff = 1:100
                 dataset(:,8) = datasample(dataset(:,8),length(dataset(:,8)),'Replace',false);   
                ds = mat2dataset(dataset,'VarNames',{'Previous','Outcome','Modality','Location','Approach','Latency','Trial','FiringRate'});
                
                %             mdl{kk}= stepwiseglm(ds,'constant','upper','linear','Distribution','poisson');
               try
                mdl{kk}.(cat(2,'shuff_',num2str(iShuff))) = stepwiseglm(ds,'constant','upper','interactions','Distribution','poisson','PEnter',.01);
           catch %err
%     switch err.identifier
%         case 'stats:glmfit:IterationLimit'
      mdl{kk}.(cat(2,'shuff_',num2str(iShuff))) = [];
%       continue
% %     elseif strcmp(err.identifier, 'Warning: Weights are ill-conditioned.   Data may be badly scaled, or the link function may be inappropriate.')
% %       mdl{kk} = [];
% %       continue
%         otherwise
%       rethrow(err); %some other unexpected error. Better stop
     end
%             'stats:glmfit:BadScaling'
                end
            end
        else
            mdl{kk} = [];
        end
    else
        mdl{kk} = [];
    end
end


save(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-03-16-GLM_cueon_',num2str(iTime),'-SHUFF.mat'),'mdl')
clearvars -except iTime
    end
%%
    for iTime = -.4:.1:.5
        load(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-03-16-GLM_cueon_',num2str(iTime),'-SHUFF.mat'))

mat_files = dir('*.mat');
%      disp('initial shuff');

for iShuff = 1:100
    
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined = [];
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc = [];
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc = [];
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec = [];
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec = [];
ALL_matrix.(cat(2,'shuff_',num2str(iShuff))) = [];

count.(cat(2,'shuff_',num2str(iShuff))) =  1;
count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))) = 1;
count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) = 1;
count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))) = 1;
count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) = 1;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Cue = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Cue = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Modality = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Modality = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Location = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Location = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Outcome = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Outcome = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Approach = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Approach = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Latency = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Latency = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Trial = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Trial = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Previous = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Previous = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLoc = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLoc = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxOut = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.LocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.LocxOut = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.OutxApp = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.OutxApp = 0;

summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLocxOut = 0;
summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLocxOut = 0;

end

for kk = 1:length(dir('*.mat'))
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
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
    
    class.all.isi{kk} = diff(spk_t);
    class.all.median_isi(kk) = median(class.all.isi{kk});
    class.all.frate(kk) = FRATE.Overall.firing_rate_total;
    class.all.sorted_isi{kk} = sort(class.all.isi{kk},'descend');
    
    if isempty(mdl{kk}) == 0
    for iShuff = 1:100
    
    switch isempty(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))))
        case 0
%              disp(cat(2,'Cell #',num2str(kk),'& Shuff #',num2str(iShuff)));
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
                    block_drift.summary(kk) = 0;
                    if RANK.two.Trial > 975 || RANK.two.Trial < 26
                        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue + 1;
                            switch class.all.sorted_isi{kk}(5) < 2
                                case 1
                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Cue + 1;
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue + 1;
                                        case 0
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Cue + 1;
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue + 1;
                                    end
                                case 0
                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Cue + 1;
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue + 1;
                                        case 0
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Cue + 1;
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue + 1;
                                    end
                            end
                            
                            for ll = 2:length(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames)
                                if mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Coefficients{ll,4} < .01 %&& mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Coefficients{ll,4} ~= 0
                                    idx = [];
                                    dev_location = [];
                                    dev_raw = [];
                                    dev_percent = [];
                                    switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Modality');
                                        case 1
                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality + 1;
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                            ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,1) = 1;
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Modality{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Modality');
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Modality{kk}.Rsquared.Adjusted) * 100;
                                            idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Modality');
                                            dev_location = find(idx == 1);
                                            dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                            dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                            switch class.all.sorted_isi{kk}(5) < 2
                                                case 1
                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality + 1;
                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                        case 1
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Modality + 1;
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality + 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                        case 0
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Modality + 1;
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality + 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                    end
                                                case 0
                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality + 1;
                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                        case 1
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Modality + 1;
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality + 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                        case 0
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Modality + 1;
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality + 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                    end
                                            end
                                        case 0
                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Location');
                                                case 1
                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location + 1;
                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,2) = 1;
                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Location{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Location');
                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Location{kk}.Rsquared.Adjusted) * 100;
                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Location');
                                                    dev_location = find(idx == 1);
                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                        case 1
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location + 1;
                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                case 1
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Location + 1;
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                case 0
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Location + 1;
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                            end
                                                        case 0
                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location + 1;
                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                case 1
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Location + 1;
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                case 0
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Location + 1;
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                            end
                                                    end
                                                case 0
                                                    switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome');
                                                        case 1
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,3) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Outcome{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Outcome{kk}.Rsquared.Adjusted) * 100;
                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome');
                                                                    dev_location = find(idx == 1);
                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                        case 1
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Outcome + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Outcome + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                            end
                                                                        case 0
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Outcome + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Outcome + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                            end
                                                                    end
                                                                case 0
                                                                 switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Approach');
                                                case 1    
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Approach + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,4) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Approach{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Approach');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Approach{kk}.Rsquared.Adjusted) * 100;
                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Approach');
                                                                    dev_location = find(idx == 1);
                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                        case 1
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Approach + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Approach + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Approach + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                            end
                                                                        case 0
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Approach + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Approach + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Approach + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                            end
                                                                    end
                                                            
                                                        case 0
                                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Latency');
                                                                case 1
                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency + 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,5) = 1;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Latency{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Latency');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Latency{kk}.Rsquared.Adjusted) * 100;
                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Latency');
                                                                    dev_location = find(idx == 1);
                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                        case 1
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Latency + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Latency + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                            end
                                                                        case 0
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency + 1;
                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Latency + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Latency + 1;
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                            end
                                                                    end
                                                                case 0
                                                                    switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Trial');
                                                                        case 1
                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial + 1;
                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                            ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,6) = 1;
                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Trial{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Trial');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Trial{kk}.Rsquared.Adjusted) * 100;
                                                                            idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Trial');
                                                                            dev_location = find(idx == 1);
                                                                            dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                            dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                            switch class.all.sorted_isi{kk}(5) < 2
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial + 1;
                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                        case 1
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Trial + 1;
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial + 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                        case 0
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Trial + 1;
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial + 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                    end
                                                                                case 0
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial + 1;
                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                        case 1
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Trial + 1;
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial + 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                        case 0
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Trial + 1;
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial + 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                    end
                                                                            end
                                                                        case 0
                                                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Previous');
                                                                                case 1
                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous + 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,7) = 1;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Previous{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Previous');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Previous{kk}.Rsquared.Adjusted) * 100;
                                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Previous');
                                                                                    dev_location = find(idx == 1);
                                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                                        case 1
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous + 1;
                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                case 1
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Previous + 1;
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous + 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                case 0
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Previous + 1;
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous + 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                            end
                                                                                        case 0
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous + 1;
                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                case 1
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Previous + 1;
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous + 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                case 0
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Previous + 1;
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous + 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                            end
                                                                                    end
                                                                                case 0
                                                                                    switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Modality:Location');
                                                                                        case 1
                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc + 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                            ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,8) = 1;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLoc{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Modality:Location');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLoc{kk}.Rsquared.Adjusted) * 100;
                                                                                            idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Modality:Location');
                                                                                            dev_location = find(idx == 1);
                                                                                            dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                            dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                            switch class.all.sorted_isi{kk}(5) < 2
                                                                                                case 1
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc + 1;
                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                        case 1
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLoc + 1;
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                        case 0
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLoc + 1;
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                    end
                                                                                                case 0
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc + 1;
                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                        case 1
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLoc + 1;
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                        case 0
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLoc + 1;
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                    end
                                                                                            end
                                                                                        case 0
                                                                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Modality');
                                                                                                case 1
                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut + 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,9) = 1;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Modality');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Modality');
                                                                                                    dev_location = find(idx == 1);
                                                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                                                        case 1
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut + 1;
                                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                case 1
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxOut + 1;
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut + 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                case 0
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxOut + 1;
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut + 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                            end
                                                                                                        case 0
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut + 1;
                                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                case 1
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxOut + 1;
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut + 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                case 0
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxOut + 1;
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut + 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                            end
                                                                                                    end
                                                                                                case 0
                                                                                                    switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Location');
                                                                                                        case 1
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                            ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,10) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.LocxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Location');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.LocxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                            idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Location');
                                                                                                    dev_location = find(idx == 1);
                                                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                            switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                case 1
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut + 1;
                                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                        case 1
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.LocxOut + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                        case 0
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.LocxOut + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                    end
                                                                                                                case 0
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut + 1;
                                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                        case 1
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.LocxOut + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                        case 0
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.LocxOut + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                    end
                                                                                                            end
                                                                                                        case 0
                                                                                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Approach');
                                                                                                        case 1
                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).All.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.OutxApp + 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),11) = 1;
                                                                                                            ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,11) = 1;
                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.OutxApp{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Approach');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,11) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.OutxApp{kk}.Rsquared.Adjusted) * 100;
                                                                                                            idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Approach');
                                                                                                    dev_location = find(idx == 1);
                                                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),11) = dev_raw;
                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),11) = dev_percent;
                                                                                                            switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                case 1
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.OutxApp + 1;
                                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                        case 1
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.OutxApp + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),11) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),11) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),11) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),11) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,11);
                                                                                                                        case 0
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.OutxApp + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,11);
                                                                                                                    end
                                                                                                                case 0
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.OutxApp + 1;
                                                                                                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                        case 1
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.OutxApp + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,11);
                                                                                                                        case 0
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.OutxApp + 1;
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp + 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = 1;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_raw;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = dev_percent;
                                                                                                                            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),11) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,11);
                                                                                                                    end
                                                                                                            end
                                                                                                        case 0
                                                                                                            switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Modality:Location');
                                                                                                                case 1
                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut + 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),12) = 1;
                                                                                                                    ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,12) = 1;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLocxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Modality:Location');
                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,12) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLocxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                                    idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Modality:Location');
                                                                                                                    dev_location = find(idx == 1);
                                                                                                                    dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                                    dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),12) = dev_raw;
                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),12) = dev_percent;
                                                                                                                    switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                        case 1
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut + 1;
                                                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                case 1
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLocxOut + 1;
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut + 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),12) = 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),12) = dev_raw;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),12) = dev_percent;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),12) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,12);
                                                                                                                                case 0
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLocxOut + 1;
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut + 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_raw;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_percent;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,12);
                                                                                                                            end
                                                                                                                        case 0
                                                                                                                           summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut + 1;
                                                                                                                            switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                case 1
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLocxOut + 1;
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut + 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_raw;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_percent;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,12);
                                                                                                                                case 0
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLocxOut + 1;
                                                                                                                                   summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut + 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = 1;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_raw;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = dev_percent;
                                                                                                                                    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),12) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,12);
                                                                                                                            end
                                                                                                                    end
                                                                                                                case 0
                                                                                                            end
                                                                                                    end
                                                                                            end
                                                                                    end
                                                                            end
                                                                    end
                                                            end
                                                    end
                                            end
                                    end
                                end
                                    end
                                end
                            end
                            count.(cat(2,'shuff_',num2str(iShuff))) =  count.(cat(2,'shuff_',num2str(iShuff))) + 1;
                            switch class.all.sorted_isi{kk}(5) < 2
                                case 1
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))) = count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                        case 0
                                            count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))) = count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                    end
                                case 0
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) = count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                        case 0
                                            count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) = count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                    end
                            end
                        end
                    end
                case 1
                    block_drift.summary(kk) = 1;
            end
    end
    end
    end
end

save(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-03-17-GLM_cueon_',num2str(iTime),'-SHUFF.mat'),'mdl','ALL_matrix','block_drift','GLM_matrices','summary_var')
clearvars -except iTime
    end

%%
summary_table{1,1} = 'Category';
summary_table{2,1} = 'Cue-evoked';
summary_table{3,1} = 'Modality';
summary_table{4,1} = 'Location';
summary_table{5,1} = 'Outcome';
summary_table{6,1} = 'Approach';
summary_table{7,1} = 'Latency';
summary_table{8,1} = 'Trial';
summary_table{9,1} = 'Previous trial';
summary_table{10,1} = 'Modality x Location';
summary_table{11,1} = 'Modality x Outcome';
summary_table{12,1} = 'Location x Outcome';
summary_table{13,1} = 'Outcome x Approach';
summary_table{14,1} = 'Modality x Location x Outcome';

summary_table{1,2} = 'All';
summary_table{2,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue;
summary_table{3,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality;
summary_table{4,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location;
summary_table{5,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome;
summary_table{6,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Approach;
summary_table{7,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency;
summary_table{8,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial;
summary_table{9,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous;
summary_table{10,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc;
summary_table{11,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut;
summary_table{12,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut;
summary_table{13,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.OutxApp;
summary_table{14,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut;

summary_table{1,3} = 'SPN';
summary_table{2,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue;
summary_table{3,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality;
summary_table{4,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location;
summary_table{5,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome;
summary_table{6,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Approach;
summary_table{7,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency;
summary_table{8,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial;
summary_table{9,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous;
summary_table{10,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc;
summary_table{11,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut;
summary_table{12,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut;
summary_table{13,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.OutxApp;
summary_table{14,3} =summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut;

summary_table{1,4} = 'HFN';
summary_table{2,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue;
summary_table{3,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality;
summary_table{4,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location;
summary_table{5,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome;
summary_table{6,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Approach;
summary_table{7,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency;
summary_table{8,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial;
summary_table{9,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous;
summary_table{10,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc;
summary_table{11,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut;
summary_table{12,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut;
summary_table{13,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.OutxApp;
summary_table{14,4} =summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut;

summary_table{1,5} = 'Exc';
summary_table{2,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue;
summary_table{3,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality;
summary_table{4,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location;
summary_table{5,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome;
summary_table{6,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Approach;
summary_table{7,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency;
summary_table{8,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial;
summary_table{9,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous;
summary_table{10,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc;
summary_table{11,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut;
summary_table{12,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut;
summary_table{13,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.OutxApp;
summary_table{14,5} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut;

summary_table{1,6} = 'Inh';
summary_table{2,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue;
summary_table{3,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality;
summary_table{4,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location;
summary_table{5,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome;
summary_table{6,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Approach;
summary_table{7,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency;
summary_table{8,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial;
summary_table{9,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous;
summary_table{10,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc;
summary_table{11,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut;
summary_table{12,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut;
summary_table{13,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.OutxApp;
summary_table{14,6} =summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut;

%% consistency check
consistency_check{1,1} = 'Category';
consistency_check{2,1} = 'Cue-evoked';
consistency_check{3,1} = 'Modality';
consistency_check{4,1} = 'Location';
consistency_check{5,1} = 'Outcome';
consistency_check{6,1} = 'Approach';
consistency_check{7,1} = 'Latency';
consistency_check{8,1} = 'Trial';
consistency_check{9,1} = 'Previous trial';
consistency_check{10,1} = 'Modality x Location';
consistency_check{11,1} = 'Modality x Outcome';
consistency_check{12,1} = 'Location x Outcome';
consistency_check{13,1} = 'Outcome x Approach';
consistency_check{14,1} = 'Modality x Location x Outcome';

consistency_check{1,2} = 'mat_overview';
consistency_check{2,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue;
consistency_check{3,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality;
consistency_check{4,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location;
consistency_check{5,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome;
consistency_check{6,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Approach;
consistency_check{7,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency;
consistency_check{8,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial;
consistency_check{9,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous;
consistency_check{10,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc;
consistency_check{11,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut;
consistency_check{12,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut;
consistency_check{13,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.OutxApp;
consistency_check{14,2} =summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut;

consistency_check{1,3} = 'GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined';
consistency_check{2,3} = [];
consistency_check{3,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,1));
consistency_check{4,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,2));
consistency_check{5,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,3));
consistency_check{6,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,4));
consistency_check{7,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,5));
consistency_check{8,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,6));
consistency_check{9,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,7));
consistency_check{10,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,8));
consistency_check{11,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,9));
consistency_check{12,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,10));
consistency_check{13,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,11));
consistency_check{14,3} = sum(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(:,12));

%%
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined_transposed = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined';

figure
imagesc(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined_transposed);
map = [229/255,245/255,249/255
    44/255,162/255,95/255];
colormap(map);
set(gca,'YTickLabel','')

%%

SPN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec * -1;
HFN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec * -1;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.cat_combined = cat(1,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc,SPN_Dec_inverse,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc,HFN_Dec_inverse)';
% cell 84 acts weird (unusually high values)
map = [233/255,163/255,201/255
    247/255,247/255,247/255
    161/255,215/255,106/255];

rColorMap = [linspace(233/255, 255/255, 111),linspace(255/255, 161/255, 145)];
    gColorMap = [linspace(163/255, 255/255, 111),linspace(255/255, 215/255, 145)];
    bColorMap = [linspace(201/255, 255/255, 111),linspace(255/255, 106/255, 145)];
colorMap = [rColorMap; gColorMap; bColorMap]';

figure
heatmap(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.cat_combined, [], [], '%0.0f', 'Colormap',colorMap, ...
        'FontSize', 0, 'Colorbar', true);
    
%%

SPN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec * -1;
HFN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec * -1;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined = cat(1,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc,SPN_Dec_inverse,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc,HFN_Dec_inverse)';
% cell 84 acts weird (unusually high values)
map = [233/255,163/255,201/255
    247/255,247/255,247/255
    161/255,215/255,106/255];

rColorMap = [linspace(233/255, 255/255, 111),linspace(255/255, 161/255, 145)];
    gColorMap = [linspace(163/255, 255/255, 111),linspace(255/255, 215/255, 145)];
    bColorMap = [linspace(201/255, 255/255, 111),linspace(255/255, 106/255, 145)];
colorMap = [rColorMap; gColorMap; bColorMap]';

figure
heatmap(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined, [], [], '%0.0f', 'Colormap',colorMap, ...
        'FontSize', 0, 'Colorbar', true);   
    
%% find average % of variance explained for each predictor
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN == 0) = NaN;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN > 50) = NaN;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN < -50) = NaN;
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Modality = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Modality = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Location = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,2));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Location = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,2))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,2))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,2))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Outcome = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,3));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Outcome = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,3))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,3))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,3))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Approach = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,4));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Approach = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,4))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,4))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,4))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Latency = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,5));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Latency = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,5))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,5))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,5))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Trial = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,6));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Trial = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,6))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,6))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,6))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Previous = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,7));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Previous = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,7))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,7))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,7))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ModxLoc = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,8));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.ModxLoc = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,8))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,8))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,8))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ModxOut = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,9));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.ModxOut = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,9))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,9))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,9))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.LocxOut = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,10));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.LocxOut = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,10))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,10))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,10))));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.OutxApp = nanmean(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,11));
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.OutxApp = nanstd(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,11))/sqrt(numel(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,11))-sum(isnan(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,11))));

GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Summary = cat(2,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Modality,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Location,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Outcome,...
    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Approach,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Latency,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Trial,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Previous);

GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Summary = cat(2,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Modality,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Location,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Outcome,...
    GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Approach,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Latency,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Trial,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Previous);

%%
figure;
bar(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Summary)
hold on
errorbar(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.Summary,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SEM.Summary,'.')
set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
ylim([0 12]);
ylabel('Percent improvement to R-Squared')

%% plot for paper
%add rows
%add (6) spaces between SPN and HFNs
GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined = cat(1,abs(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc),abs(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec)*-1,abs(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc),abs(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec)*-1)';
for iCell = length(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined):-1:1
    for iFeature = 1:length(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined(:,iCell))
    if GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined(iFeature,iCell) > 100 || GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined(iFeature,iCell) < -100
        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined(:,iCell) = [];
        break
    end
    end
end

%%
rColorMap = [linspace(233/255, 255/255, 183),linspace(255/255, 161/255, 73)];
    gColorMap = [linspace(163/255, 255/255, 183),linspace(255/255, 215/255, 73)];
    bColorMap = [linspace(201/255, 255/255, 183),linspace(255/255, 106/255, 73)];
colorMap = [rColorMap; gColorMap; bColorMap]';

figure
subplot(3,12,[2 3 4 5 6 14 15 16 17 18])
heatmap(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined(1:10,:),  'RowLabels', {'Cue identity','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue identity x location','Cue identity x outcome','Cue location x outcome'},... 
'%0.0f', 'Colormap',colorMap, ...
        'FontSize', 0); %, 'Colorbar', true);   
% set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
xlabel('Unit number');
% xlabel('MSNs increasing      MSNs decreasing          FSIs increasing        FSIs decreasing')
title('GLM matrix for cue-modulated units');  
set(gca,'FontSize',18);
hold on; plot(108,1:.001:10,'.k');
subplot(3,12,[8 9 10 11 12 20 21 22 23 24])
violin(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1:7));
set(gca,'XTickLabel',{'Identity','Location','Outcome','Approach','Trial length','Trial number','Previous'},'XTickLabelRotation',45)
% ylim([0 12]);
ylabel('Percent improvement to R-Squared')
title('Summary of changes to model fit')
set(gca,'FontSize',18);

%%

% SPN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec * -1;
% HFN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec * -1;
% GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.cat_combined = cat(1,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc,SPN_Dec_inverse,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc,HFN_Dec_inverse)';
% % cell 84 acts weird (unusually high values)
% map = [233/255,163/255,201/255
%     247/255,247/255,247/255
%     161/255,215/255,106/255];
% 
% rColorMap = [linspace(233/255, 255/255, 111),linspace(255/255, 161/255, 145)];
%     gColorMap = [linspace(163/255, 255/255, 111),linspace(255/255, 215/255, 145)];
%     bColorMap = [linspace(201/255, 255/255, 111),linspace(255/255, 106/255, 145)];
% colorMap = [rColorMap; gColorMap; bColorMap]';
% 
% figure
% heatmap(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.cat_combined(1:10,:), 'RowLabels', {'Cue modality','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue modality x location','Cue modality x outcome','Cue location x outcome'},... 
% '%0.0f', 'Colormap',colorMap,'FontSize', 0, 'Colorbar', true);
% xlabel('Unit number');
% % xlabel('MSNs increasing                                            MSNs decreasing                        FSIs increasing   FSIs decreasing')
% title('GLM matrix for cue-modulated units');        
%%
% 
% SPN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec * -1;
% HFN_Dec_inverse = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec * -1;
% GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined = cat(1,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc,SPN_Dec_inverse,GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc,HFN_Dec_inverse)';
% % cell 84 acts weird (unusually high values)
% map = [233/255,163/255,201/255
%     247/255,247/255,247/255
%     161/255,215/255,106/255];
% 
% rColorMap = [linspace(233/255, 255/255, 190),linspace(255/255, 161/255, 66)];
%     gColorMap = [linspace(163/255, 255/255, 190),linspace(255/255, 215/255, 66)];
%     bColorMap = [linspace(201/255, 255/255, 190),linspace(255/255, 106/255, 66)];
% colorMap = [rColorMap; gColorMap; bColorMap]';
% 
% figure
% heatmap(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.cat_combined,  'RowLabels', {'Cue modality','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous trial','Cue modality x location','Cue modality x outcome','Cue location x outcome'},... 
% '%0.0f', 'Colormap',colorMap, ...
%         'FontSize', 0, 'Colorbar', true);   
% set(gca,'XTickLabel',{'Mod','Loc','Out','App','Lat','Trial','Prev'})
% xlabel('Unit number');
% % xlabel('MSNs increasing      MSNs decreasing          FSIs increasing        FSIs decreasing')
%     
% 
% %% violin plot
% figure;
% violin(GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN.ALL_NaN(:,1:7));
% set(gca,'XTickLabel',{'Cue modality','Cue location','Cue outcome','Approach','Trial length','Trial number','Previous'})
% % ylim([0 12]);
% ylabel('Percent improvement to R-Squared')