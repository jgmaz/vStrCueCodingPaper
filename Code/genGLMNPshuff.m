function GLMnpshuff = genGLMNPshuff(directory,destination,num_Shuffs)
% function GLMnpshuff = genGLMNPshuff(directory,destination,num_Shuffs)

warning('error', 'stats:glmfit:IterationLimit')
rng('shuffle')
%% time window to analyze
for iEpoch = 1:2
    for iTime = -.5:.1:.5
        
        time_window_start = iTime; %starting time window for analysis, 0 = time zero
        
        switch iEpoch %which epoch to look at, 1 = np delay, 2 = outcome receipt
            case 1
                time_window_end = time_window_start + .5;
                epoch_start = 0;
                Epoch = 'NP';
            case 2
                time_window_start = time_window_start + 1;
                time_window_end = time_window_start + .5;
                epoch_start = 1;
                Epoch = 'outcome';
        end
        
        %%
        mat_files = dir('*.mat');
        for kk = 1:length(dir('*.mat'))
            load(mat_files(kk).name);
            mat_overview.fname{kk} = mat_files(kk).name;
            disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
            
            if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                
                clear dataset ds %mdl
                trials_count = 1;
                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                for jj = 1:b1_length
                    if metadata.TrialInfo{1,1}.nosepoke_length(jj) > epoch_start
                        if metadata.TrialInfo{1,1}.nosepoke_length(jj) < time_window_end
                            nosepoke_to_click(trials_count) = metadata.TrialInfo{1,1}.nosepoke_length(jj);
                        else
                            nosepoke_to_click(trials_count) = time_window_end;
                        end
                        end_time(trials_count) = metadata.dataPoint.Nosepokes(trials_count)/1000 + nosepoke_to_click(trials_count);
                        start_time(trials_count) = metadata.dataPoint.Nosepokes(trials_count)/1000 + time_window_start;
                        
                        % now, count.(cat(2,'shuff_',num2str(iShuff))) spikes between start and end
                        these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                        
                        % convert to firing rate and store
                        firing_rate_np(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
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
                        dataset(trials_count,5) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                        dataset(trials_count,6) = jj;
                        dataset(trials_count,7) = firing_rate_np(trials_count); %FR
                        
                        trials_count = trials_count + 1;
                    end
                end
                
                b2_start = b1_length + 1;
                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                total = b1_length + b2_length;
                
                for jj = 1:b2_length
                    if metadata.TrialInfo{1,2}.nosepoke_length(jj) > epoch_start
                        if metadata.TrialInfo{1,2}.nosepoke_length(jj) < time_window_end
                            nosepoke_to_click(trials_count) = metadata.TrialInfo{1,2}.nosepoke_length(jj);
                        else
                            nosepoke_to_click(trials_count) = time_window_end;
                        end
                        end_time(trials_count) = metadata.dataPoint.Nosepokes(trials_count)/1000 + nosepoke_to_click(trials_count);
                        start_time(trials_count) = metadata.dataPoint.Nosepokes(trials_count)/1000 + time_window_start;
                        
                        % now, count.(cat(2,'shuff_',num2str(iShuff))) spikes between start and end
                        these_spk = spk_t(spk_t > start_time(trials_count) & spk_t < end_time(trials_count));
                        
                        % convert to firing rate and store
                        firing_rate_np(trials_count) = length(these_spk) / (end_time(trials_count) - start_time(trials_count));
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
                        dataset(trials_count,5) = metadata.TrialInfo{1,2}.summary(jj,15);
                        dataset(trials_count,6) = b1_length+jj;
                        dataset(trials_count,7) = firing_rate_np(trials_count);
                        
                        trials_count = trials_count + 1;
                    end
                end
                
                if sum(dataset(:,7)) < 30
                    mdl{kk} = [];
                else
                    %%
                    for iShuff = 1:num_Shuffs
                        dataset(:,7) = datasample(dataset(:,7),length(dataset(:,7)),'Replace',false);
                        ds = mat2dataset(dataset,'VarNames',{'Previous','Outcome','Modality','Location','Latency','Trial','FiringRate'});
                        try
                            mdl{kk}.(cat(2,'shuff_',num2str(iShuff))) = stepwiseglm(ds,'constant','upper','interactions','Distribution','poisson','PEnter',.01);
                        catch %err
                            mdl{kk}.(cat(2,'shuff_',num2str(iShuff))) = [];
                        end
                    end
                end
            else
                mdl{kk} = [];
            end
        end
        
% save(cat(2,destination,'GLM_',Epoch,'_SHUFF_',num2str(iTime),'-intermediate_file.mat'),'mdl')
        
%         clearvars -except iEpoch iTime
%     end
% end
%%
% Epoch = {'NP' 'outcome'};
% for iEpoch = 2
%     for iTime = -.5:.1:.5
%         load(cat(2,destination,'GLM_',Epoch,'_SHUFF_',num2str(iTime),'-intermediate_file.mat'))
        
        mat_files = dir('*.mat');
        
        for iShuff = 1:num_Shuffs
            
            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined = [];
            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc = [];
            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc = [];
            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec = [];
            GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec = [];
            ALL_matrix.(cat(2,'shuff_',num2str(iShuff))) = [];
            
            count.(cat(2,'shuff_',num2str(iShuff))) = 1;
            count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))) = 1;
            count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))) = 1;
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
            class.all.sorted_isi{kk} = sort(class.all.isi{kk},'descend');
            
            if isempty(mdl{kk}) == 0
                for iShuff = 1:num_Shuffs
                    
                    switch isempty(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))))
                        case 0
                            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                                case 0
                                    if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Cue + 1;
                                        switch class.all.sorted_isi{kk}(5) < 2
                                            case 1
                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Cue + 1;
                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                    case 1
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Cue + 1;
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue + 1;
                                                    case 0
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Cue + 1;
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue + 1;
                                                end
                                            case 0
                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Cue + 1;
                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                    case 1
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Cue + 1;
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Cue + 1;
                                                    case 0
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Cue + 1;
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Cue + 1;
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
                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Modality + 1;
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
                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Modality + 1;
                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                    case 1
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Modality + 1;
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality + 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                                    case 0
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Modality + 1;
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality + 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                                end
                                                            case 0
                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Modality + 1;
                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                    case 1
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Modality + 1;
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Modality + 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                                    case 0
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Modality + 1;
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Modality + 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = 1;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_raw;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = dev_percent;
                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),1) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,1);
                                                                end
                                                        end
                                                    case 0
                                                        switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Location');
                                                            case 1
                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Location + 1;
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
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Location + 1;
                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                            case 1
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Location + 1;
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location + 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                            case 0
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Location + 1;
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location + 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                        end
                                                                    case 0
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Location + 1;
                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                            case 1
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Location + 1;
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Location + 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                            case 0
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Location + 1;
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Location + 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_raw;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = dev_percent;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),2) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,2);
                                                                        end
                                                                end
                                                            case 0
                                                                switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome');
                                                                    case 1
                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Outcome + 1;
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
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Outcome + 1;
                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                    case 1
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Outcome + 1;
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome + 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                    case 0
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Outcome + 1;
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome + 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                end
                                                                            case 0
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Outcome + 1;
                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                    case 1
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Outcome + 1;
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Outcome + 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                    case 0
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Outcome + 1;
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Outcome + 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_raw;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = dev_percent;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),3) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,3);
                                                                                end
                                                                        end
                                                                    case 0
                                                                    case 0
                                                                        switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Latency');
                                                                            case 1
                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Latency + 1;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,4) = 1;
                                                                                idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Latency');
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Latency{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Latency');
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Latency{kk}.Rsquared.Adjusted) * 100;
                                                                                dev_location = find(idx == 1);
                                                                                dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                switch class.all.sorted_isi{kk}(5) < 2
                                                                                    case 1
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Latency + 1;
                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                            case 1
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Latency + 1;
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency + 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                            case 0
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Latency + 1;
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency + 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                        end
                                                                                    case 0
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Latency + 1;
                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                            case 1
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Latency + 1;
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Latency + 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                            case 0
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Latency + 1;
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Latency + 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_raw;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = dev_percent;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),4) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,4);
                                                                                        end
                                                                                end
                                                                            case 0
                                                                                switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Trial');
                                                                                    case 1
                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Trial + 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                        ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,5) = 1;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Trial{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Trial');
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Trial{kk}.Rsquared.Adjusted) * 100;
                                                                                        idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Trial');
                                                                                        dev_location = find(idx == 1);
                                                                                        dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                        dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                        switch class.all.sorted_isi{kk}(5) < 2
                                                                                            case 1
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Trial + 1;
                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                    case 1
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Trial + 1;
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial + 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                                    case 0
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Trial + 1;
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial + 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                                end
                                                                                            case 0
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Trial + 1;
                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                    case 1
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Trial + 1;
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Trial + 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                                    case 0
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Trial + 1;
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Trial + 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_raw;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = dev_percent;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),5) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,5);
                                                                                                end
                                                                                        end
                                                                                    case 0
                                                                                        switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Previous');
                                                                                            case 1
                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.Previous + 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                                ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,6) = 1;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Previous{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Previous');
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.Previous{kk}.Rsquared.Adjusted) * 100;
                                                                                                idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Previous');
                                                                                                dev_location = find(idx == 1);
                                                                                                dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                                switch class.all.sorted_isi{kk}(5) < 2
                                                                                                    case 1
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.Previous + 1;
                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                            case 1
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.Previous + 1;
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous + 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                                            case 0
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.Previous + 1;
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous + 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                                        end
                                                                                                    case 0
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.Previous + 1;
                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                            case 1
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.Previous + 1;
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.Previous + 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                                            case 0
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.Previous + 1;
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.Previous + 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_raw;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = dev_percent;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),6) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,6);
                                                                                                        end
                                                                                                end
                                                                                            case 0
                                                                                                switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Modality:Location');
                                                                                                    case 1
                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLoc + 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                        ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,7) = 1;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLoc{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Modality:Location');
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLoc{kk}.Rsquared.Adjusted) * 100;
                                                                                                        idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Modality:Location');
                                                                                                        dev_location = find(idx == 1);
                                                                                                        dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                        dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                        switch class.all.sorted_isi{kk}(5) < 2
                                                                                                            case 1
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLoc + 1;
                                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                    case 1
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLoc + 1;
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc + 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                                    case 0
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLoc + 1;
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc + 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                                end
                                                                                                            case 0
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLoc + 1;
                                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                    case 1
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLoc + 1;
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLoc + 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                                    case 0
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLoc + 1;
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLoc + 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_raw;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = dev_percent;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),7) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,7);
                                                                                                                end
                                                                                                        end
                                                                                                    case 0
                                                                                                        switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Modality');
                                                                                                            case 1
                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxOut + 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                                ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,8) = 1;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Modality');
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                                idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Modality');
                                                                                                                dev_location = find(idx == 1);
                                                                                                                dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                                dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                                switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                    case 1
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxOut + 1;
                                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                            case 1
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxOut + 1;
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut + 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                                            case 0
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxOut + 1;
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut + 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                                        end
                                                                                                                    case 0
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxOut + 1;
                                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                            case 1
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxOut + 1;
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxOut + 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                                            case 0
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxOut + 1;
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxOut + 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_raw;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = dev_percent;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),8) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,8);
                                                                                                                        end
                                                                                                                end
                                                                                                            case 0
                                                                                                                switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Location');
                                                                                                                    case 1
                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.LocxOut + 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                        ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,9) = 1;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.LocxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Location');
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.LocxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                                        idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Location');
                                                                                                                        dev_location = find(idx == 1);
                                                                                                                        dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                                        dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                        switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                            case 1
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.LocxOut + 1;
                                                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                    case 1
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.LocxOut + 1;
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut + 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                                    case 0
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.LocxOut + 1;
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut + 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                                end
                                                                                                                            case 0
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.LocxOut + 1;
                                                                                                                                switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                    case 1
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.LocxOut + 1;
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.LocxOut + 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                                    case 0
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.LocxOut + 1;
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.LocxOut + 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = 1;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_raw;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = dev_percent;
                                                                                                                                        GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),9) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,9);
                                                                                                                                end
                                                                                                                        end
                                                                                                                    case 0
                                                                                                                        switch strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).CoefficientNames(ll),'Outcome:Modality:Location');
                                                                                                                            case 1
                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.ModxLocxOut + 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                                ALL_matrix.(cat(2,'shuff_',num2str(iShuff)))(kk,10) = 1;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLocxOut{kk} = removeTerms(mdl{kk}.(cat(2,'shuff_',num2str(iShuff))),'Outcome:Modality:Location');
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10) = (mdl{kk}.(cat(2,'shuff_',num2str(iShuff))).Rsquared.Adjusted - GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Comparison.ModxLocxOut{kk}.Rsquared.Adjusted) * 100;
                                                                                                                                idx = strcmp(mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.TermName,'Outcome:Modality:Location');
                                                                                                                                dev_location = find(idx == 1);
                                                                                                                                dev_raw = mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1) - mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location);
                                                                                                                                dev_percent = (dev_raw / mdl{1,kk}.(cat(2,'shuff_',num2str(iShuff))).Steps.History.Deviance(dev_location-1)) * 100;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.combined(count.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                                switch class.all.sorted_isi{kk}(5) < 2
                                                                                                                                    case 1
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN.ModxLocxOut + 1;
                                                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                            case 1
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Inc.ModxLocxOut + 1;
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut + 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Inc(count_HFN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                                            case 0
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).HFN_Dec.ModxLocxOut + 1;
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut + 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.HFN_Dec(count_HFN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                                        end
                                                                                                                                    case 0
                                                                                                                                        summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN.ModxLocxOut + 1;
                                                                                                                                        switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                                                                                                                            case 1
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Inc.ModxLocxOut + 1;
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Inc.ModxLocxOut + 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Inc(count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
                                                                                                                                            case 0
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).SPN_Dec.ModxLocxOut + 1;
                                                                                                                                                summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut = summary_var.(cat(2,'shuff_',num2str(iShuff))).Dec.ModxLocxOut + 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).count.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = 1;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_raw.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_raw;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).dev_percent.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = dev_percent;
                                                                                                                                                GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.SPN_Dec(count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))),10) = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL(kk,10);
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
                                        count.(cat(2,'shuff_',num2str(iShuff))) = count.(cat(2,'shuff_',num2str(iShuff))) + 1;
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
                                                        count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))) = count_SPN_Inc.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                                    case 0
                                                        count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) = count_SPN_Dec.(cat(2,'shuff_',num2str(iShuff))) + 1;
                                                end
                                        end
                                    end
                            end
                    end
                end
            end
        end
        
        save(cat(2,destination,'GLM_',Epoch,'_SHUFF_',num2str(iTime),'.mat'),'mdl','ALL_matrix','block_drift','GLM_matrices','summary_var')
        clearvars -except iTime destination directory num_Shuffs iEpoch
    end
end

end