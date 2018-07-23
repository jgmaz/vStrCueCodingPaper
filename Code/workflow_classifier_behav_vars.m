%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Incorporating location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get data into matrix for LDA


mat_files = dir('*.mat');
t_count = 0;

for iTime = -.5:.1:.5
    t_count = t_count + 1;
    count = 1;
    cell_count = 1;
    
    LDA_raw_input{t_count}(1:500,1:133) = NaN;
    
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
                        
                        
%                         switch sesh.block_order %modality
%                             case 1
%                                 starting_b1 = 0;
%                                 starting_b2 = 120;
%                             case 2
%                                 starting_b1 = 120;
%                                 starting_b2 = 0;
%                         end
starting_b1 = 0;
starting_b2 = 0;
                        
                        switch new_v_old
                            case 0
                                %% old rats (R053,R056,R057)
                                trials_count = 1;
                                short_app_count = 301;
                                long_app_count = 201;
                                short_skip_count = 101;
                                long_skip_count = 1;
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
                                        
                                        switch metadata.TrialInfo_block1.approached(jj)
                                            case 0
                                                if metadata.TrialInfo_block1.trial_length_analysis(jj) < mean(metadata.TrialInfo_block1.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_skip_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        short_skip_count = short_skip_count + 1;
                                                else
                                                        LDA_raw_input{t_count}(long_skip_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        long_skip_count = long_skip_count + 1;
                                                end
                                            case 1
                                                if metadata.TrialInfo_block1.trial_length_analysis(jj) < mean(metadata.TrialInfo_block1.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_app_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        short_app_count = short_app_count + 1;
                                                else
                                                        LDA_raw_input{t_count}(long_app_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        long_app_count = long_app_count + 1;
                                                end
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                
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
                                        
                                        switch metadata.TrialInfo_block2.approached(jj)
                                            case 0
                                                if metadata.TrialInfo_block2.trial_length_analysis(jj) < mean(metadata.TrialInfo_block2.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_skip_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        short_skip_count = short_skip_count + 1;
                                                        else
                                                        LDA_raw_input{t_count}(long_skip_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        long_skip_count = long_skip_count + 1;
                                                end
                                            case 1
                                                if metadata.TrialInfo_block2.trial_length_analysis(jj) < mean(metadata.TrialInfo_block2.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_app_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        short_app_count = short_app_count + 1;
                                                else
                                                        LDA_raw_input{t_count}(long_app_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        long_app_count = long_app_count + 1;
                                                end
                                        end
                                        
                                        trials_count = trials_count + 1;

                                end
                                
                            case 1
                                %% new rats (R060)
                                trials_count = 1;
                                short_app_count = 301;
                                long_app_count = 201;
                                short_skip_count = 101;
                                long_skip_count = 1;
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
                                    
                                    switch metadata.TrialInfo{1,1}.approached(jj)
                                        case 0
                                            if metadata.TrialInfo{1,1}.trial_length_analysis(jj) < mean(metadata.TrialInfo{1,1}.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_skip_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        short_skip_count = short_skip_count + 1;
                                            else
                                                        LDA_raw_input{t_count}(long_skip_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        long_skip_count = long_skip_count + 1;
                                                end
                                        case 1
                                            if metadata.TrialInfo{1,1}.trial_length_analysis(jj) < mean(metadata.TrialInfo{1,1}.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_app_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        short_app_count = short_app_count + 1;
                                            else
                                                        LDA_raw_input{t_count}(long_app_count+starting_b1,cell_count) = firing_rate(trials_count);
                                                        long_app_count = long_app_count + 1;
                                                end
                                    end
                                    
                                    trials_count = trials_count + 1;
                                end
                                
                                b2_start = b1_length + 1;
                                b2_length = length(FRATE.Cue.Trial_firing_rate_block2);
                                total = b1_length + b2_length;
                                
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
                                    
                                    switch metadata.TrialInfo{1,2}.approached(jj)
                                        case 0
                                            if metadata.TrialInfo{1,2}.trial_length_analysis(jj) < mean(metadata.TrialInfo{1,2}.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_skip_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        short_skip_count = short_skip_count + 1;
                                                        else
                                                        LDA_raw_input{t_count}(long_skip_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        long_skip_count = long_skip_count + 1;
                                                end
                                        case 1
                                            if metadata.TrialInfo{1,2}.trial_length_analysis(jj) < mean(metadata.TrialInfo{1,2}.trial_length_analysis)
                                                        LDA_raw_input{t_count}(short_app_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        short_app_count = short_app_count + 1;
                                            else
                                                        LDA_raw_input{t_count}(long_app_count+starting_b2,cell_count) = firing_rate(trials_count);
                                                        long_app_count = long_app_count + 1;
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
mdl_identifier = {'Latency' 'Approach' 'AppxLat'};

for iLabel = 1:100    
    Labels.Latency{iLabel,1} = 'Long';
    Labels.Approach{iLabel,1} = 'Skip';
    Labels.AppxLat{iLabel,1} = 'LongSkip';
end
for iLabel = 101:200
    Labels.Latency{iLabel,1} = 'Short';
    Labels.Approach{iLabel,1} = 'Skip';
    Labels.AppxLat{iLabel,1} = 'ShortSkip';
end
for iLabel = 201:300
    Labels.Latency{iLabel,1} = 'Long';
    Labels.Approach{iLabel,1} = 'App';
    Labels.AppxLat{iLabel,1} = 'LongApp';
end
for iLabel = 301:400
    Labels.Latency{iLabel,1} = 'Short';
    Labels.Approach{iLabel,1} = 'App';
    Labels.AppxLat{iLabel,1} = 'ShortApp';
end


%% remove cells with 0 variance for each label
for iTime = 1:length(LDA_raw_input) %note, skipped 4 and 5
    for iMdl = 1:length(mdl_identifier)
        disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))
        Label_loop = unique(Labels.(mdl_identifier{iMdl}));
        Label_exclude = [];
        for iLabel = 1:length(Label_loop)
           Label_trials = strcmp(Labels.(mdl_identifier{iMdl}),Label_loop{iLabel});
               Label_variance = nanvar(LDA_raw_input{iTime}(Label_trials,:));
               Label_exclude = cat(2,Label_exclude,find(Label_variance < .6)); %note used .6 here as a cut-off instead of .15 to get LDA to work consistently
        end
        LDA_fixed_input.(mdl_identifier{iMdl}){iTime} = LDA_raw_input{iTime};
        LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(:,Label_exclude) = [];       
    end
end
       
%% Cross-validation
rng('shuffle')

for iTime = 1:length(LDA_raw_input) %note, skipped 4 and 5
    for iMdl = 1:length(mdl_identifier)
        disp(cat(2,'Mdl ',num2str(iMdl),' (',num2str(iTime),')'))        
    for iteration = 1:25
        LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(1:400,:),Labels.(mdl_identifier{iMdl}));
%     switch iMdl
%         case 1
%         LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(1:400,[1:5 7:17 19:24 26:end]),Labels.(mdl_identifier{iMdl}));
%         case 2
%              LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(1:400,[2:15 17:18 20:end]),Labels.(mdl_identifier{iMdl}));
%         case 3
%              LDA_Mdl = fitcdiscr(LDA_fixed_input.(mdl_identifier{iMdl}){iTime}(1:400,[1:5 7:16 19:end]),Labels.(mdl_identifier{iMdl}));
%     end
        LDA_cv = crossval(LDA_Mdl);
    LDA_Error.(mdl_identifier{iMdl})(iteration,iTime) = kfoldLoss(LDA_cv);
    end   
    end
end

% [1:91 93:96 98:124 126:133]
%% Convert MCR to CR
for iMdl = 1:length(mdl_identifier)
Class_accuracy.(mdl_identifier{iMdl}) = 1 - LDA_Error.(mdl_identifier{iMdl});
end

%% Plot LDA
figure
for iMdl = 1:length(mdl_identifier)
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