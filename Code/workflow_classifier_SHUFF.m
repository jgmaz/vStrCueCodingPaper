%% get data into matrix for LDA


mat_files = dir('*.mat');
t_count = 0;

for iTime = -.5:.1:.5
    cell_count = 1;
    t_count = t_count + 1;
    count = 1;
    LDA_all_trials{t_count}(1:200,1:133) = NaN;
    
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
                                starting_b2 = 100;
                            case 2
                                starting_b1 = 100;
                                starting_b2 = 0;
                        end
                        
                        switch new_v_old
                            case 0
                                %% old rats (R053,R056,R057)
                                trials_count = 1;
                                rew_count = 1;
                                unrew_count = 51;
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    if jj > 100
                                        break
                                    end
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
                                        
                                        switch metadata.TrialInfo_block1.rewarded(jj)
                                            case 0
                                                LDA_all_trials{t_count}(unrew_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                unrew_count = unrew_count + 1;
                                            case 1
                                                LDA_all_trials{t_count}(rew_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                rew_count = rew_count + 1;
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                rew_count = 1;
                                unrew_count = 51;
                                
                                for jj = 1:b2_length
                                   if jj > 100
                                       break
                                   end
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
                                        
                                        switch metadata.TrialInfo_block2.rewarded(jj)
                                            case 0
                                                LDA_all_trials{t_count}(unrew_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                unrew_count = unrew_count + 1;
                                            case 1
                                                LDA_all_trials{t_count}(rew_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                rew_count = rew_count + 1;
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                            case 1
                                %% new rats (R060)
                                trials_count = 1;
                                rew_count = 1;
                                unrew_count = 51;
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    if jj > 100
                                        break
                                    end
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
                                    
                                    switch metadata.TrialInfo{1,1}.rewarded(jj)
                                        case 0
                                            LDA_all_trials{t_count}(unrew_count+starting_b1,cell_count) = firing_rate(trials_count);
                                            unrew_count = unrew_count + 1;
                                        case 1
                                            LDA_all_trials{t_count}(rew_count+starting_b1,cell_count) = firing_rate(trials_count);
                                            rew_count = rew_count + 1;
                                    end
                                    
                                    trials_count = trials_count + 1;
                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                rew_count = 1;
                                unrew_count = 51;
                                
                                for jj = 1:b2_length
                                    if jj > 100
                                        break
                                    end
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
                                    
                                    switch metadata.TrialInfo{1,2}.rewarded(jj)
                                        case 0
                                            LDA_all_trials{t_count}(unrew_count+starting_b2,cell_count) = firing_rate(trials_count);
                                            unrew_count = unrew_count + 1;
                                        case 1
                                            LDA_all_trials{t_count}(rew_count+starting_b2,cell_count) = firing_rate(trials_count);
                                            rew_count = rew_count + 1;
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

%% Cross-val, how to get plot?
for iLabel = 1:25
Labels.ModxOut{iLabel,1} = 'L+';
Labels.Modality{iLabel,1} = 'L';
Labels.Outcome{iLabel,1} = '+';
Labels.Drift{iLabel,1} = '1st';
end
for iLabel = 26:50
Labels.ModxOut{iLabel,1} = 'L+';
Labels.Modality{iLabel,1} = 'L';
Labels.Outcome{iLabel,1} = '+';
Labels.Drift{iLabel,1} = '2nd';
end
for iLabel = 51:75
Labels.ModxOut{iLabel,1} = 'L-';
Labels.Modality{iLabel,1} = 'L';
Labels.Outcome{iLabel,1} = '-';
Labels.Drift{iLabel,1} = '1st';
end
for iLabel = 76:100
Labels.ModxOut{iLabel,1} = 'L-';
Labels.Modality{iLabel,1} = 'L';
Labels.Outcome{iLabel,1} = '-';
Labels.Drift{iLabel,1} = '2nd';
end
for iLabel = 101:125
Labels.ModxOut{iLabel,1} = 'S+';
Labels.Modality{iLabel,1} = 'S';
Labels.Outcome{iLabel,1} = '+';
Labels.Drift{iLabel,1} = '1st';
end
for iLabel = 126:150
Labels.ModxOut{iLabel,1} = 'S+';
Labels.Modality{iLabel,1} = 'S';
Labels.Outcome{iLabel,1} = '+';
Labels.Drift{iLabel,1} = '2nd';
end
for iLabel = 151:175
Labels.ModxOut{iLabel,1} = 'S-';
Labels.Modality{iLabel,1} = 'S';
Labels.Outcome{iLabel,1} = '-';
Labels.Drift{iLabel,1} = '2nd';
end
for iLabel = 176:200
Labels.ModxOut{iLabel,1} = 'S-';
Labels.Modality{iLabel,1} = 'S';
Labels.Outcome{iLabel,1} = '-';
Labels.Drift{iLabel,1} = '2nd';
end

%%
mdl_identifier = {'Modality' 'Outcome' 'ModxOut' 'Drift' 'Drift_L' 'Drift_S'};
for iTime = 1:length(LDA_all_trials) %note, skipped 4 and 5
    for iMdl = 5%6 %1:length(mdl_identifier)
        disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
    for iteration = 1:100
    LDA_Mdl = fitcdiscr(LDA_all_trials{iTime}(1:100,[1:91 93:96 98:124 126:133]),Labels.(mdl_identifier{iMdl}));
    LDA_cv = crossval(LDA_Mdl);
    LDA_Error.(mdl_identifier{iMdl})(iteration,iTime) = kfoldLoss(LDA_cv);
    end   
    end
end

%%
for iMdl = 1:length(mdl_identifier)
Class_accuracy.(mdl_identifier{iMdl}) = 1 - LDA_Error.(mdl_identifier{iMdl});
end

%%
figure
for iMdl = 1:length(mdl_identifier)
subplot(2,2,iMdl)
violin(Class_accuracy.(mdl_identifier{iMdl}));
set(gca,'XTickLabel',{'','','','','','','','','','',''})
ylim([0 1])
hold on;
switch iMdl
    case 1
plot(-.5:.025:.5,.5,'.','color','black');
ylabel('Classification accuracy')
    case 2
plot(-.5:.025:.5,.5,'.','color','black');
    case 3
plot(-.5:.025:.5,.25,'.','color','black');
ylabel('Classification accuracy')
xlabel('Pseudoensemble size')
    case 4
plot(-.5:.025:.5,.5,'.','color','black');

xlabel('Pseudoensemble size')
end
title(mdl_identifier{iMdl})
set(gca,'FontSize',18);
end

%% LDA
for iLabel = 1:49
LDA_label{iLabel,1} = 'L+';
end
for iLabel = 50:98
LDA_label{iLabel,1} = 'L-';
end
for iLabel = 99:147
LDA_label{iLabel,1} = 'S+';
end
for iLabel = 148:196
LDA_label{iLabel,1} = 'S-';
end


LDA_train = cat(1,LDA_all_trials{t_count}(2:50,:),LDA_all_trials{t_count}(52:100,:),LDA_all_trials{t_count}(102:150,:),LDA_all_trials{t_count}(152:200,:));
LDA_test = cat(1,LDA_all_trials{t_count}(1,:),LDA_all_trials{t_count}(51,:),LDA_all_trials{t_count}(101,:),LDA_all_trials{t_count}(151,:));
Mdl = fitcdiscr(LDA_train,LDA_label);
LDA_error = resubLoss(Mdl);
LDA_predict = predict(Mdl,LDA_test);

%% Cross-val, how to get plot?
for iLabel = 1:50
LDA_label_all{iLabel,1} = 'L+';
end
for iLabel = 51:100
LDA_label_all{iLabel,1} = 'L-';
end
for iLabel = 101:150
LDA_label_all{iLabel,1} = 'S+';
end
for iLabel = 151:200
LDA_label_all{iLabel,1} = 'S-';
end

Mdl_all = fitcdiscr(LDA_all_trials,LDA_label_all);
cv_Mdl = crossval(Mdl_all);
cv_Error = kfoldLoss(cv_Mdl);

%%
LDA_remove = LDA_all_trials{t_count}(:,[1:91 93:124 126:133]); 

rng('shuffle');
count = 0;
for iSize = 2:6:133
    count = count + 1;
    for iteration = 1:100
        disp(cat(2,'iteration: ',num2str(iteration),' (Size: ',num2str(iSize),')'))
        use_cells = randsample(131,iSize);
    LDA_subset = LDA_remove(:,use_cells);
    LDA_Mdl = fitcdiscr(LDA_subset,LDA_label_all);
    LDA_cv = crossval(LDA_Mdl,'leaveout','on');
    LDA_Error(iteration,count) = kfoldLoss(LDA_cv);
    end    
end

%%
Class_accuracy = 1 - LDA_Error;
Class_avg = mean(Class_accuracy);

%%
figure;
violin(Class_accuracy');
set(gca,'XTickLabel',{'8','20','32','44','56','68','80','92','104','116','128'})
ylim([0 1])
hold on;
plot(-.5:.025:.5,.25,'.','color','black');
ylabel('Classification accuracy')
xlabel('Pseudoensemble size')
title('Classification accuracy as a function of Pseudoensemble size')
set(gca,'FontSize',18);

%% Cross-val, out or mod
for iLabel = 1:50
LDA_label_out{iLabel,1} = '+';
end
for iLabel = 51:100
LDA_label_out{iLabel,1} = '-';
end
for iLabel = 101:150
LDA_label_out{iLabel,1} = '+';
end
for iLabel = 151:200
LDA_label_out{iLabel,1} = '-';
end

    for iteration = 1:100
        disp(cat(2,'iteration: ',num2str(iteration),' (Size: ',num2str(iSize),')'))
%         use_cells = randsample(131,iSize);
%     LDA_subset = LDA_remove(:,use_cells);
    LDA_Mdl = fitcdiscr(LDA_all_trials,LDA_label_out);
    LDA_cv = crossval(LDA_Mdl);
    LDA_Error(iteration) = kfoldLoss(LDA_cv);
    end
    
%% Cross-val, out or mod
for iLabel = 1:100
LDA_label_mod{iLabel,1} = 'L';
end

for iLabel = 101:200
LDA_label_mod{iLabel,1} = 'S';
end

    for iteration = 1:100
        disp(cat(2,'iteration: ',num2str(iteration),' (Size: ',num2str(iSize),')'))
%         use_cells = randsample(131,iSize);
%     LDA_subset = LDA_remove(:,use_cells);
    LDA_Mdl = fitcdiscr(LDA_all_trials,LDA_label_mod);
    LDA_cv = crossval(LDA_Mdl);
    LDA_Error(iteration) = kfoldLoss(LDA_cv);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incorporating location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get data into matrix for LDA
rng('shuffle')

mat_files = dir('*.mat');
t_count = 0;
num_Shuffs = 100;

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
                                
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    %                                     if jj > 100
                                    %                                         break
                                    %                                     end
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
                                    
                                    trials_count = trials_count + 1;
                                end
                                
                                
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                
                                
                                for jj = 1:b2_length
                                    %                                    if jj > 100
                                    %                                        break
                                    %                                    end
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
                                        
                                        switch metadata.TrialInfo_block1.rewarded(jj)
                                            case 0
                                                switch metadata.TrialInfo_block1.photosensorID(jj)
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
                                                switch metadata.TrialInfo_block1.photosensorID(jj)
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
                                        switch metadata.TrialInfo_block2.rewarded(jj)
                                            case 0
                                                switch metadata.TrialInfo_block2.photosensorID(jj)
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
                                                switch metadata.TrialInfo_block2.photosensorID(jj)
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
                                
                            case 1
                                %% new rats (R060)
                                trials_count = 1;
                                
                                b1_length = length(FRATE.Cue.Trial_firing_rate_block1);
                                for jj = 1:b1_length
                                    %                                     if jj > 100
                                    %                                         break
                                    %                                     end
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
                                    
                                    trials_count = trials_count + 1;
                                end
                                
                                
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                
                                
                                for jj = 1:b2_length
                                    %                                    if jj > 100
                                    %                                        break
                                    %                                    end
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

%%
for iTime = 1:length(LDA_raw_input) 
    for iShuff = 1:100
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
    for iShuff = 1:100
        for iMdl = 1:length(mdl_identifier)
            disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
            for iteration = 1:100
                LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}.(cat(2,'s',num2str(iShuff)))(1:240,:),Labels.(mdl_identifier{iMdl}));
                LDA_cv = crossval(LDA_Mdl);
                LDA_Error.(mdl_identifier{iMdl}).(cat(2,'s',num2str(iShuff)))(iteration,iTime) = kfoldLoss(LDA_cv);
            end
        end
    end    
end

%% 
for iShuff = 1:100
    for iMdl = 1:length(mdl_identifier)
        LDA_Error_SHUFF.(mdl_identifier{iMdl})(iShuff,:) = nanmean(LDA_Error.(mdl_identifier{iMdl}).(cat(2,'s',num2str(iShuff))));
    end
end

%%
for iMdl = 1:length(mdl_identifier)
Class_accuracy_SHUFF.(mdl_identifier{iMdl}) = 1 - LDA_Error_SHUFF.(mdl_identifier{iMdl});
end

%%
figure
for iMdl = 1:length(mdl_identifier)
subplot(3,3,iMdl)
% violin(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
hold on
% violin(Class_accuracy.(mdl_identifier{iMdl}));
% set(gca,'XTickLabel',{'','','','','','','','','','',''})
shadedErrorBar(-.5:.1:.5,mean(Class_accuracy.(mdl_identifier{iMdl})),std(Class_accuracy.(mdl_identifier{iMdl}))/sqrt(length(Class_accuracy_SHUFF.(mdl_identifier{iMdl}))),'-b',1);
ylim([0 1])

shadedErrorBar(-.5:.1:.5,mean(Class_accuracy_SHUFF.(mdl_identifier{iMdl})),std(Class_accuracy_SHUFF.(mdl_identifier{iMdl}))/sqrt(length(Class_accuracy_SHUFF.(mdl_identifier{iMdl}))),'-r',1);
switch iMdl
    case 1
plot(-.5:.025:.5,.5,'.','color','black');
    case 2
plot(-.5:.025:.5,.25,'.','color','black');
    case 3
plot(-.5:.025:.5,.5,'.','color','black');
    case 4
plot(-.5:.025:.5,.25,'.','color','black');
ylabel('Classification accuracy')
    case 5
plot(-.5:.025:.5,.125,'.','color','black');
    case 6
plot(-.5:.025:.5,.125,'.','color','black');
    case 7
plot(-.5:.025:.5,.0625,'.','color','black');
xlabel('Time')
end
title(mdl_identifier{iMdl})
set(gca,'FontSize',18);
for iTime = 1:length(Table.(mdl_identifier{iMdl}).Zscore)
    if Table.(mdl_identifier{iMdl}).Zscore(iTime) > 1.96
        plot((-.6+(.1*iTime)),.9,'*k')
    end
end
end

%% Z-scores away from shuffle (value - shuff mean / shuff std)
mdl_identifier = {'Modality','Location','Outcome','ModxOut','ModxLoc','LocxOut','ModxLocxOut'};

for iMdl = 1:length(mdl_identifier);
    Table.(mdl_identifier{iMdl}).Data = mean(Class_accuracy.(mdl_identifier{iMdl}));
    Table.(mdl_identifier{iMdl}).ShuffMEAN = mean(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    Table.(mdl_identifier{iMdl}).ShuffSTD = std(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    for iTime = 1:length(Table.(mdl_identifier{iMdl}).Data)
         Table.(mdl_identifier{iMdl}).Zscore(iTime) = (Table.(mdl_identifier{iMdl}).Data(iTime) - Table.(mdl_identifier{iMdl}).ShuffMEAN(iTime)) / Table.(mdl_identifier{iMdl}).ShuffSTD(iTime);
            if Table.(mdl_identifier{iMdl}).Zscore(iTime) > 1.96
                Table.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 1;
            else
                Table.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 0;
            end
    end
end

%%
figure
colors = {'b' 'r' 'y'};
for iEpoch = 1%:length(Epochs)
    subplot(1,3,iEpoch)
    hold on
for iMdl = 1:length(mdl_identifier);
    plot(-.5:.1:.5,Table.(mdl_identifier{iMdl}).Zscore,'color',colors{iMdl})
end
title(Epochs{iEpoch})
plot(-.5:.05:.5,1.96,'.','color','black');
ylabel('Zscore')
xlabel('Time')
end

%%
for iMdl = 1:length(mdl_identifier);
    Class_accuracy.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy.(mdl_identifier{iMdl}));
    Class_accuracy.SEM.(mdl_identifier{iMdl}) = std(Class_accuracy.(mdl_identifier{iMdl}))/sqrt(length(Class_accuracy.(mdl_identifier{iMdl})));
    Class_accuracy_SHUFF.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    Class_accuracy_SHUFF.SEM.(mdl_identifier{iMdl}) = std(Class_accuracy_SHUFF.(mdl_identifier{iMdl}))/sqrt(length(Class_accuracy_SHUFF.(mdl_identifier{iMdl})));
end
    
%% cue on w SHUFF
figure;
% Mdls = {[1 4 5 7] [2 5 6 7] [3 4 6 7]};
Mdls = {[1 2 3] [4 5 6] [7]};
% colors = {{'b' 'r' 'y' 'g'}}; % {'r' 'b' 'y' 'g'} {'y' 'b' 'r' 'g'}};
colors = {'b' 'r' 'y' 'g'};
start_pt = [.05 .05 .03];

for iMdl = 1:length(Mdls)

subplot(2,3,iMdl)
hold on

for iPlot = 1:length(Mdls{iMdl})
shadedErrorBar(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy.SEM.(mdl_identifier{Mdls{iMdl}(iPlot)}),colors{iPlot},1);
shadedErrorBar(-.5:.1:.5,Class_accuracy_SHUFF.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy_SHUFF.SEM.(mdl_identifier{Mdls{iMdl}(iPlot)}),strcat('--',colors{iPlot}),1);
for iTime = 1:length(Table.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode)
    if Table.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[start_pt(iMdl)-(iPlot*.01) start_pt(iMdl)-(iPlot*.01)], '-k', 'LineWidth',2,'color',colors{iPlot})
    end
end
end
plot(-.05,0.01:.01:1,'.k'); plot(-.45,0.01:.01:1,'.k');
% legend({'identity' 'location' 'outcome'}); 
switch iMdl
    case 1        
title('Cue features');
ylabel('Classification rate')
    case 2
        title('Multiple cue features')
    case 3
        title('All cue features')
end

ylim([0 1]); xlabel('LDA start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);

end

%%
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.SEM.Outcome,'--y',1);


shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.SEM.Modality,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.SEM.Outcome,'--y',1);
for iTime = 1:length(Table.cueon.Outcome.Zscore_recode)
    if Table.cueon.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.cueon.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.cueon.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','y')
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);

subplot(2,3,2)
plot(-.5:.1:.5,GLM_window.cueon.prop(5:6,:))
hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Approach,GLM_window_SHUFF.cueon.prop.SEM.Approach,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Latency,GLM_window_SHUFF.cueon.prop.SEM.Latency,'--r',1);
for iTime = 1:length(Table.cueon.Outcome.Zscore_recode)
    if Table.cueon.Approach.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
    end
    if Table.cueon.Latency.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
legend({'approach' 'latency'});

subplot(2,3,3)
plot(-.5:.1:.5,GLM_window.cueon.prop(7:8,:))
hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Trial,GLM_window_SHUFF.cueon.prop.SEM.Trial,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Previous,GLM_window_SHUFF.cueon.prop.SEM.Previous,'--r',1);
for iTime = 1:length(Table.cueon.Outcome.Zscore_recode)
    if Table.cueon.Trial.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
    end
    if Table.cueon.Previous.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
legend({'trial number' 'previous trial'});

subplot(2,3,4)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window_SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Location,GLM_window_SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window_SHUFF.cueon.Rsquared.SEM.Outcome,'--y',1);
ylim([1 8]); 
  y = ylabel('Percent improvement to R-Squared');
set(y, 'position', get(y,'position')+[-.0005,0,0]); 
xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
subplot(2,3,5)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(4,:),GLM_window.cueon.Rsquared.SEM(4,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(5,:),GLM_window.cueon.Rsquared.SEM(5,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Approach,GLM_window_SHUFF.cueon.Rsquared.SEM.Approach,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Latency,GLM_window_SHUFF.cueon.Rsquared.SEM.Latency,'--r',1);
ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
subplot(2,3,6)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(6,:),GLM_window.cueon.Rsquared.SEM(6,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(7,:),GLM_window.cueon.Rsquared.SEM(7,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Trial,GLM_window_SHUFF.cueon.Rsquared.SEM.Trial,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Previous,GLM_window_SHUFF.cueon.Rsquared.SEM.Previous,'--r',1);
ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);