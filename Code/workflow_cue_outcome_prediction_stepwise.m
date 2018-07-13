ds = mat2dataset(cat(2,dataset(:,3),dataset(:,3)),'VarNames',{'Outcome','Approach'});

%             mdl{t_count,iSesh}= stepwiseglm(ds,'constant','upper','linear','Distribution','poisson');
%                try
mdl{t_count,iSesh}= stepwiseglm(ds,'constant','upper','interactions','Distribution','binomial','PEnter',.01);
%            catch %err


tbl = table(dataset(:,2),dataset(:,3),'VariableNames',{'oneback','Outcome'});
tbl.Outcome = categorical(tbl.Outcome);
tbl.oneback = categorical(tbl.oneback);
glm_test = fitglm(tbl,'Link','logit','Distribution','binomial');
glm_linear = fitglm(tbl);
figure
hist(glm_test.Residuals.Raw,100)

out = glm_linear.Fitted;
out = categorical(out.Response > 0.5);
data = categorical(dataset(:,3) > 0.5);
pred = sum(out == data)/length(dataset(:,3)); % count number of correct predictions

tbl = table(dataset(:,1),dataset(:,2),dataset(:,3),'VariableNames',{'twoback','oneback','Outcome'});
tbl.Outcome = categorical(tbl.Outcome);
tbl.oneback = categorical(tbl.oneback);
tbl.twoback = categorical(tbl.twoback);
glm_test = fitglm(tbl,'Link','logit','Distribution','binomial');
glm_linear2 = fitglm(tbl);
glm_linear_int = fitglm(tbl,'interactions');

tbl = table(dataset(:,9),dataset(:,3),'VariableNames',{'frate','Outcome'});
tbl.Outcome = categorical(tbl.Outcome);
glm_test = fitglm(tbl,'Link','logit','Distribution','binomial');
glm_linear = fitglm(tbl);

tbl = table(dataset(:,9),dataset(:,1),dataset(:,2),dataset(:,3),'VariableNames',{'frate','twoback','oneback','Outcome'});
tbl.Outcome = categorical(tbl.Outcome);
tbl.oneback = categorical(tbl.oneback);
tbl.twoback = categorical(tbl.twoback);
glm_test = fitglm(tbl,'Link','logit','Distribution','binomial');
glm_linear = fitglm(tbl);

glm_linear2.Deviance

%% model comparisons

%% time window to analyze
warning('error', 'stats:glmfit:IterationLimit'); % turn GLM warning into error so that you can use try catch to skip to next iteration of loop.

% modelspec = {'Outcome ~ 1 + Previous','Outcome ~ 1 + PrevTwo','Outcome ~ 1 + PrevApp' ...
%     'Outcome ~ 1 + RRate','Outcome ~ 1 + Previous + PrevTwo','Outcome ~ 1 + Previous + PrevApp', ...
%     'Outcome ~ 1 + Previous + RRate','Outcome ~ 1'};
% mdl_identifier = {'OneBack','TwoBack','App','RRate','OneTwoBack','OneApp','OneRRate','Int'};
% modelspec = {'Outcome ~ 1 + Previous','Outcome ~ 1 + FiringRate','Outcome ~ 1 + Previous + FiringRate','Outcome ~ 1'};
% mdl_identifier = {'OneBack','FRate','OneFRate','Int'};
modelspec = {'Outcome ~ 1 + Previous','Outcome ~ 1 + PrevTwo','Outcome ~ 1 + PrevApp' ...
    'Outcome ~ 1 + RRate','Outcome ~ 1 + FiringRate','Outcome ~ 1 + Previous + FiringRate', ...
    'Outcome ~ 1 + PrevApp + FiringRate','Outcome ~ 1 + RRate + FiringRate','Outcome ~ 1'};
mdl_identifier = {'OneBack','TwoBack','App','RRate','FRate','OneFRate','AppFRate','RRateFRate','Int'};

for iMdl = 1:length(mdl_identifier)
    count{iMdl} = 1;
    count_bdrift{iMdl} = 1;
    count_cuemod{iMdl} = 1;
end

t_count = 0;
for iTime = -.5:.1:.5
    t_count = t_count + 1;
    
    for iMdl = 1:length(mdl_identifier)
        count{iMdl} = 1;
        count_bdrift{iMdl} = 1;
        count_cuemod{iMdl} = 1;
    end
    
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    mat_files = dir('*.mat');
    for iSesh = 1:length(dir('*.mat'))
        load(mat_files(iSesh).name);
        mat_overview.fname{iSesh} = mat_files(iSesh).name;
        disp(cat(2,num2str(iSesh),'/',num2str(length(dir('*.mat')))));
        new_v_old = strcmp(mat_overview.fname{iSesh}(1:4),'R060');
        
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo_block1.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                        
                    else
                        dataset(trials_count,1) = metadata.TrialInfo_block1.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo_block1.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                    switch sesh.block_order %modality
                        case 1
                            dataset(trials_count,4) = 1;
                        case 2
                            dataset(trials_count,4) = 2;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo_block1.summary(jj,5); %arm
                    dataset(trials_count,6) = metadata.TrialInfo_block1.approached(jj); %behav - app
                    dataset(trials_count,7) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                    dataset(trials_count,8) = jj;
                    dataset(trials_count,9) = firing_rate(trials_count); %FR
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo_block2.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo_block2.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo_block2.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo_block2.rewarded(jj);
                    switch sesh.block_order
                        case 1
                            dataset(trials_count,4) = 2;
                        case 2
                            dataset(trials_count,4) = 1;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo_block2.summary(jj,5);
                    dataset(trials_count,6) = metadata.TrialInfo_block2.approached(jj);
                    dataset(trials_count,7) = metadata.TrialInfo_block2.summary(jj,15);
                    dataset(trials_count,8) = b1_length+jj;
                    dataset(trials_count,9) = firing_rate(trials_count);
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo{1,1}.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo{1,1}.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo{1,1}.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                    switch sesh.block_order %modality
                        case 1
                            dataset(trials_count,4) = 1;
                        case 2
                            dataset(trials_count,4) = 2;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                    dataset(trials_count,6) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                    dataset(trials_count,7) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                    dataset(trials_count,8) = jj;
                    dataset(trials_count,9) = firing_rate(trials_count); %FR
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo{1,2}.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo{1,2}.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo{1,2}.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo{1,2}.rewarded(jj);
                    switch sesh.block_order
                        case 1
                            dataset(trials_count,4) = 2;
                        case 2
                            dataset(trials_count,4) = 1;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo{1,2}.summary(jj,5);
                    dataset(trials_count,6) = metadata.TrialInfo{1,2}.approached(jj);
                    dataset(trials_count,7) = metadata.TrialInfo{1,2}.summary(jj,15);
                    dataset(trials_count,8) = b1_length+jj;
                    dataset(trials_count,9) = firing_rate(trials_count);
                    
                    trials_count = trials_count + 1;
                    
                end
                
        end
        
        %         ds = mat2dataset(dataset,'VarNames',{'PrevTwo','Previous','Outcome','Modality','Location','Approach','Latency','Trial','FiringRate','PrevApp','RRate'});
        dataset2 = cat(2,dataset(:,2),dataset(:,9),dataset(:,3));
        ds = mat2dataset(dataset2,'VarNames',{'Previous','FRate','Outcome'});
        ds.Outcome = categorical(ds.Outcome);
        ds.Previous = categorical(ds.Previous);
        
        block_drift.block1_length(iSesh) = length(FRATE.Cue.Trial_firing_rate_block1);
        block_drift.block1_half(iSesh) = round(block_drift.block1_length(iSesh) / 2);
        block_drift.b1_1st_avg(iSesh) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(iSesh)));
        block_drift.b1_2nd_avg(iSesh) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(iSesh)+1:end));
        block_drift.MWU_b1(iSesh) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(iSesh)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(iSesh)+1:end));
        
        block_drift.block2_length(iSesh) = length(FRATE.Cue.Trial_firing_rate_block2);
        block_drift.block2_half(iSesh) = round(block_drift.block2_length(iSesh) / 2);
        block_drift.b2_1st_avg(iSesh) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(iSesh)));
        block_drift.b2_2nd_avg(iSesh) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(iSesh)+1:end));
        block_drift.MWU_b2(iSesh) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(iSesh)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(iSesh)+1:end));
        
        switch block_drift.MWU_b1(iSesh) < .01 || block_drift.MWU_b2(iSesh) < .01
            case 0
                bdrift = 1;
            case 1
                bdrift = 0;
        end
        
        switch TESTS.WSR.Task.Trial_b4_vs_Trial < .01
            case 0
                cuemod = 0;
            case 1
                cuemod = 1;
        end
        iMdl = 1;
        %         for iMdl = 1:length(modelspec)
        try
            mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9},'Link','logit','Distribution','binomial');
            MdlComp{t_count}.ALL.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Deviance;
            MdlComp{t_count}.ALL.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp{t_count}.ALL.Deviance(iSesh,1)-MdlComp{t_count}.ALL.Deviance(iSesh,iMdl),1);
            out = mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Fitted;
            out = categorical(out.Response > 0.5);
            data = categorical(dataset(:,3) > 0.5);
            MdlComp{t_count}.ALL.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
            MdlComp{t_count}.ALL.predRatio(iSesh,iMdl) = MdlComp{t_count}.ALL.predicted(iSesh,iMdl) / MdlComp{t_count}.ALL.predicted(iSesh,1);
            for iCoeff = 1:length(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames)
                if strcmp(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FRate') == 1
                    if mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                        MdlComp{t_count}.ALL.FRate(iSesh,iMdl) = 1;
                        MdlComp{t_count}.ALL.Count = count{iMdl};
                        MdlComp{t_count}.ALL.FRate_cells{iMdl}(count{iMdl}) = MdlComp{t_count}.ALL.predicted(iSesh,iMdl);
                        count{iMdl} = count{iMdl} + 1;
                    end
                end
            end
            
            if bdrift == 1
                mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9},'Link','logit','Distribution','binomial');
                MdlComp{t_count}.bdrift.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Deviance;
                MdlComp{t_count}.bdrift.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp{t_count}.bdrift.Deviance(iSesh,1)-MdlComp{t_count}.bdrift.Deviance(iSesh,iMdl),1);
                out = mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Fitted;
                out = categorical(out.Response > 0.5);
                data = categorical(dataset(:,3) > 0.5);
                MdlComp{t_count}.bdrift.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
                MdlComp{t_count}.bdrift.predRatio(iSesh,iMdl) = MdlComp{t_count}.bdrift.predicted(iSesh,iMdl) / MdlComp{t_count}.bdrift.predicted(iSesh,1);
                for iCoeff = 1:length(mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).CoefficientNames)
                    if strcmp(mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FRate') == 1
                        if mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                            MdlComp{t_count}.bdrift.FRate(iSesh,iMdl) = 1;
                            MdlComp{t_count}.bdrift.Count = count_bdrift{iMdl};
                            MdlComp{t_count}.bdrift.FRate_cells{iMdl}(count_bdrift{iMdl}) = MdlComp{t_count}.bdrift.predicted(iSesh,iMdl);
                            count_bdrift{iMdl} = count_bdrift{iMdl} + 1;
                        else
                            MdlComp{t_count}.bdrift.FRate(iSesh,iMdl) = 0;
                        end
                    end
                end
                
                if cuemod == 1
                    mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9},'Link','logit','Distribution','binomial');
                    MdlComp{t_count}.cuemod.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Deviance;
                    MdlComp{t_count}.cuemod.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp{t_count}.cuemod.Deviance(iSesh,1)-MdlComp{t_count}.cuemod.Deviance(iSesh,iMdl),1);
                    out = mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Fitted;
                    out = categorical(out.Response > 0.5);
                    data = categorical(dataset(:,3) > 0.5);
                    MdlComp{t_count}.cuemod.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
                    MdlComp{t_count}.cuemod.predRatio(iSesh,iMdl) = MdlComp{t_count}.cuemod.predicted(iSesh,iMdl) / MdlComp{t_count}.cuemod.predicted(iSesh,1);
                    for iCoeff = 1:length(mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).CoefficientNames)
                        if strcmp(mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FRate') == 1
                            if mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                                MdlComp{t_count}.cuemod.FRate(iSesh,iMdl) = 1;
                                MdlComp{t_count}.cuemod.Count = count_cuemod{iMdl};
                                MdlComp{t_count}.cuemod.FRate_cells{iMdl}(count_cuemod{iMdl}) = MdlComp{t_count}.cuemod.predicted(iSesh,iMdl);
                                count_cuemod{iMdl} = count_cuemod{iMdl} + 1;
                            else
                                MdlComp{t_count}.cuemod.FRate(iSesh,iMdl) = 0;
                            end
                        end
                    end
                end
            end
            
            
        catch
            mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}) = [];
        end
        %         end
        
        %         MdlComp{t_count}.ALL.Proportions(iSesh,1) = sum(dataset(:,2) == 0 & dataset(:,3) == 0);
        %         MdlComp{t_count}.ALL.Proportions(iSesh,2) = sum(dataset(:,2) == 1 & dataset(:,3) == 0);
        %         MdlComp{t_count}.ALL.Proportions(iSesh,3) = sum(dataset(:,2) == 0 & dataset(:,3) == 1);
        %         MdlComp{t_count}.ALL.Proportions(iSesh,4) = sum(dataset(:,2) == 1 & dataset(:,3) == 1);
        %         MdlComp{t_count}.ALL.CondP(iSesh,1) = MdlComp{t_count}.ALL.Proportions(iSesh,1) / (MdlComp{t_count}.ALL.Proportions(iSesh,1) + MdlComp{t_count}.ALL.Proportions(iSesh,2));
        %         MdlComp{t_count}.ALL.CondP(iSesh,2) = MdlComp{t_count}.ALL.Proportions(iSesh,2) / (MdlComp{t_count}.ALL.Proportions(iSesh,1) + MdlComp{t_count}.ALL.Proportions(iSesh,2));
        %         MdlComp{t_count}.ALL.CondP(iSesh,3) = MdlComp{t_count}.ALL.Proportions(iSesh,3) / (MdlComp{t_count}.ALL.Proportions(iSesh,3) + MdlComp{t_count}.ALL.Proportions(iSesh,4));
        %         MdlComp{t_count}.ALL.CondP(iSesh,4) = MdlComp{t_count}.ALL.Proportions(iSesh,4) / (MdlComp{t_count}.ALL.Proportions(iSesh,3) + MdlComp{t_count}.ALL.Proportions(iSesh,4));
    end
    
    clearvars -except iTime mdl MdlComp t_count modelspec mdl_identifier
end

% ds.Outcome = categorical(ds.Outcome);
% ds.Previous = categorical(ds.Previous);
% modelspec = 'Outcome ~ 1 + Previous';
% glm = fitglm(ds,modelspec);
%
% p = 1-chi2cdf(glm_linear.Deviance-glm_linear2.Deviance,1)

%%
MdlComp{t_count}.ALL.FRate_count = 0;

for iSesh = 1:length(mdl)
    for iCoeff = 1:length(mdl{t_count,iSesh}.ALL.OneBack.CoefficientNames)
        if strcmp(mdl{t_count,iSesh}.ALL.OneBack.CoefficientNames(iCoeff),'FRate') == 1
            if mdl{t_count,iSesh}.ALL.OneBack.Coefficients{iCoeff,4} < .01
                MdlComp{t_count}.ALL.FRate(iSesh) = 1;
                MdlComp{t_count}.ALL.FRate_count = MdlComp{t_count}.ALL.FRate_count + 1;
            else
                MdlComp{t_count}.ALL.FRate(iSesh) = 0;
            end
        end
    end
end

%%
MdlComp{t_count}.bdrift.FRate_count = 0;

for iSesh = 1:length(mdl)
    if isfield(mdl{t_count,iSesh},'bdrift') == 1
        for iCoeff = 1:length(mdl{t_count,iSesh}.bdrift.OneBack.CoefficientNames)
            if strcmp(mdl{t_count,iSesh}.bdrift.OneBack.CoefficientNames(iCoeff),'FRate') == 1
                if mdl{t_count,iSesh}.bdrift.OneBack.Coefficients{iCoeff,4} < .01
                    MdlComp{t_count}.bdrift.FRate(iSesh) = 1;
                    MdlComp{t_count}.bdrift.FRate_count = MdlComp{t_count}.bdrift.FRate_count + 1;
                else
                    MdlComp{t_count}.bdrift.FRate(iSesh) = 0;
                end
            end
        end
    end
end

%%
for iMdl = 1:length(MdlComp{t_count}.ALL.FRate(1,:))
    MdlComp{t_count}.ALL.SummaryRAW.FRate(iMdl) = sum(MdlComp{t_count}.ALL.FRate(:,iMdl) == 1);
    MdlComp{t_count}.bdrift.SummaryRAW.FRate(iMdl) = sum(MdlComp{t_count}.bdrift.FRate(:,iMdl) == 1);
    MdlComp{t_count}.cuemod.SummaryRAW.FRate(iMdl) = sum(MdlComp{t_count}.cuemod.FRate(:,iMdl) == 1);
    
    for iMdl2 = 1:length(MdlComp{t_count}.ALL.FRate(1,:))
        MdlComp{t_count}.ALL.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp{t_count}.ALL.FRate(:,iMdl) == 1 & MdlComp{t_count}.ALL.FRate(:,iMdl2) == 1);
        MdlComp{t_count}.ALL.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp{t_count}.ALL.FRate(:,iMdl) == 1 & MdlComp{t_count}.ALL.FRate(:,iMdl2) == 0);
        
        MdlComp{t_count}.bdrift.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp{t_count}.bdrift.FRate(:,iMdl) == 1 & MdlComp{t_count}.bdrift.FRate(:,iMdl2) == 1);
        MdlComp{t_count}.bdrift.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp{t_count}.bdrift.FRate(:,iMdl) == 1 & MdlComp{t_count}.bdrift.FRate(:,iMdl2) == 0);
        
        MdlComp{t_count}.cuemod.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp{t_count}.cuemod.FRate(:,iMdl) == 1 & MdlComp{t_count}.cuemod.FRate(:,iMdl2) == 1);
        MdlComp{t_count}.cuemod.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp{t_count}.cuemod.FRate(:,iMdl) == 1 & MdlComp{t_count}.cuemod.FRate(:,iMdl2) == 0);
    end
end


%%
count_beta = 1;
for iCell = 1:length(MdlComp{1,1}.ALL.FRate)
    if MdlComp{1,1}.ALL.FRate(iCell) == 1
        for iCoeff = 1:length(mdl{1,iCell}.ALL.OneBack.CoefficientNames)
            if strcmp(mdl{1,iCell}.ALL.OneBack.CoefficientNames(iCoeff),'FRate') == 1
                beta_FRate(count_beta) = mdl{1,iCell}.ALL.OneBack.Coefficients{iCoeff,1}
                count_beta = count_beta + 1;
            end
        end
    end
end
%%
count_all = 1;
count_bd = 1;
count_cm = 1;

for iCell = 1:442%length(mdl)
    if MdlComp{t_count}.ALL.FRate(iCell,5) == 1
        sigF.ALL(count_all,:) = MdlComp{t_count}.ALL.predicted(iCell,:);
        count_all = count_all + 1;
    end
    if MdlComp{t_count}.bdrift.FRate(iCell,5) == 1
        sigF.bd(count_bd,:) = MdlComp{t_count}.bdrift.predicted(iCell,:);
        count_bd = count_bd + 1;
    end
    if MdlComp{t_count}.cuemod.FRate(iCell,5) == 1
        sigF.cm(count_cm,:) = MdlComp{t_count}.cuemod.predicted(iCell,:);
        count_cm = count_cm + 1;
    end
end

% MdlComp{t_count}.ALL.SummaryRAW.Combined = sum(MdlComp{t_count}.ALL.FRate(:,3) == 1);
% MdlComp{t_count}.ALL.SummaryRAW.FRateOnly = sum(MdlComp{t_count}.ALL.FRate(:,2) == 1 & MdlComp{t_count}.ALL.FRate(:,3) == 0);
% MdlComp{t_count}.ALL.SummaryRAW.CombinedOnly = sum(MdlComp{t_count}.ALL.FRate(:,3) == 1 & MdlComp{t_count}.ALL.FRate(:,2) == 0);
% MdlComp{t_count}.ALL.SummaryRAW.Both = sum(MdlComp{t_count}.ALL.FRate(:,2) == 1 & MdlComp{t_count}.ALL.FRate(:,3) == 1);
% MdlComp{t_count}.ALL.SummaryPROP.FRate = sum(MdlComp{t_count}.ALL.FRate(:,2) == 1) /133;
% MdlComp{t_count}.ALL.SummaryPROP.Combined = sum(MdlComp{t_count}.ALL.FRate(:,3) == 1) /133;
% MdlComp{t_count}.ALL.SummaryPROP.FRateOnly = sum(MdlComp{t_count}.ALL.FRate(:,2) == 1 & MdlComp{t_count}.ALL.FRate(:,3) == 0) /133;
% MdlComp{t_count}.ALL.SummaryPROP.CombinedOnly = sum(MdlComp{t_count}.ALL.FRate(:,3) == 1 & MdlComp{t_count}.ALL.FRate(:,2) == 0) /133;
% MdlComp{t_count}.ALL.SummaryPROP.Both = sum(MdlComp{t_count}.ALL.FRate(:,2) == 1 & MdlComp{t_count}.ALL.FRate(:,3) == 1) /133;
%%

% count = 1;
% for iCell = 1:length(MdlComp{t_count}.ALL.FRate)
%     if MdlComp{t_count}.ALL.FRate(iCell,2) == 1
%         FRate_cells(count) = MdlComp{t_count}.ALL.predicted(iCell,2);
%         count = count + 1;
%     end
% end
%
% count = 1;
% for iCell = 1:length(MdlComp{t_count}.ALL.FRate)
%     if MdlComp{t_count}.ALL.FRate(iCell,3) == 1
%         FRate_cells2(count) = MdlComp{t_count}.ALL.predicted(iCell,3);
%         count = count + 1;
%     end
% end
%%
mat_files = dir('*.mat');
count = 1;
iTime = -.5;
time_window_start = iTime; %starting time window for analysis, 0 = time zero
time_window_end = time_window_start + .5;
epoch_start = 0;
for iCell = 1:length(mdl)
    mdl_present = isempty(mdl{1,iCell}.ALL.OneBack);
    if mdl_present == 0
        FR_present = 0;
        Coefficients = [];
        Coefficients(1) = mdl{1,iCell}.ALL.OneBack.Coefficients{1,1};
        for iCoeff = 2:length(mdl{1,iCell}.ALL.OneBack.CoefficientNames)
            if strcmp(mdl{1,iCell}.ALL.OneBack.CoefficientNames(iCoeff),'Previous_1') == 1
                Coefficients(2) = mdl{1,iCell}.ALL.OneBack.Coefficients{iCoeff,1};
            end
            if strcmp(mdl{1,iCell}.ALL.OneBack.CoefficientNames(iCoeff),'FRate') == 1
                if mdl{1,iCell}.ALL.OneBack.Coefficients{iCoeff,4} < .01
                    FR_present = 1;
                    Coefficients(3) = mdl{1,iCell}.ALL.OneBack.Coefficients{iCoeff,1};
                end
            end
        end
        
        if FR_present == 1;
            
            load(mat_files(iCell).name);
            mat_overview.fname{iCell} = mat_files(iCell).name;
            disp(cat(2,num2str(iCell),'/',num2str(length(dir('*.mat')))));
            new_v_old = strcmp(mat_overview.fname{iCell}(1:4),'R060');
            
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo_block1.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                            
                        else
                            dataset(trials_count,1) = metadata.TrialInfo_block1.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo_block1.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                        switch sesh.block_order %modality
                            case 1
                                dataset(trials_count,4) = 1;
                            case 2
                                dataset(trials_count,4) = 2;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo_block1.summary(jj,5); %arm
                        dataset(trials_count,6) = metadata.TrialInfo_block1.approached(jj); %behav - app
                        dataset(trials_count,7) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                        dataset(trials_count,8) = jj;
                        dataset(trials_count,9) = firing_rate(trials_count); %FR
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo_block2.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo_block2.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo_block2.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo_block2.rewarded(jj);
                        switch sesh.block_order
                            case 1
                                dataset(trials_count,4) = 2;
                            case 2
                                dataset(trials_count,4) = 1;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo_block2.summary(jj,5);
                        dataset(trials_count,6) = metadata.TrialInfo_block2.approached(jj);
                        dataset(trials_count,7) = metadata.TrialInfo_block2.summary(jj,15);
                        dataset(trials_count,8) = b1_length+jj;
                        dataset(trials_count,9) = firing_rate(trials_count);
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo{1,1}.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo{1,1}.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo{1,1}.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                        switch sesh.block_order %modality
                            case 1
                                dataset(trials_count,4) = 1;
                            case 2
                                dataset(trials_count,4) = 2;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                        dataset(trials_count,6) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                        dataset(trials_count,7) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                        dataset(trials_count,8) = jj;
                        dataset(trials_count,9) = firing_rate(trials_count); %FR
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo{1,2}.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo{1,2}.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo{1,2}.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo{1,2}.rewarded(jj);
                        switch sesh.block_order
                            case 1
                                dataset(trials_count,4) = 2;
                            case 2
                                dataset(trials_count,4) = 1;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo{1,2}.summary(jj,5);
                        dataset(trials_count,6) = metadata.TrialInfo{1,2}.approached(jj);
                        dataset(trials_count,7) = metadata.TrialInfo{1,2}.summary(jj,15);
                        dataset(trials_count,8) = b1_length+jj;
                        dataset(trials_count,9) = firing_rate(trials_count);
                        
                        trials_count = trials_count + 1;
                        
                    end
                    
            end
            
            %         ds = mat2dataset(dataset,'VarNames',{'PrevTwo','Previous','Outcome','Modality','Location','Approach','Latency','Trial','FiringRate','PrevApp','RRate'});
            dataset2 = cat(2,dataset(:,2),dataset(:,9),dataset(:,3));
            
            Cue_Pred.coeff(count,:) = Coefficients;
            %     Cue_Pred.pred{count}(:,1) = dataset(:,3);
            %     Cue_Pred.pred{count}(:,2) = Cue_Pred.coeff(count,1) + (dataset(:,2) * Cue_Pred.coeff(count,2));
            %     Cue_Pred.pred{count}(:,3) = Cue_Pred.coeff(count,1) + (dataset(:,2) * Cue_Pred.coeff(count,2)) + (dataset(:,9) * Cue_Pred.coeff(count,3));
            
            Cue_Pred.pred{count}(:,1) = dataset(:,3);
            for iTrial = 1:length(dataset(:,2))
                Cue_Pred.pred{count}(iTrial,2) = exp(Cue_Pred.coeff(count,1) + (dataset(iTrial,2) * Cue_Pred.coeff(count,2))) / (1 + exp(Cue_Pred.coeff(count,1) + (dataset(iTrial,2) * Cue_Pred.coeff(count,2))));
                Cue_Pred.pred{count}(iTrial,3) = exp(Cue_Pred.coeff(count,1) + (dataset(iTrial,2) * Cue_Pred.coeff(count,2)) + (dataset(iTrial,9) * Cue_Pred.coeff(count,3))) / (1 + exp(Cue_Pred.coeff(count,1) + (dataset(iTrial,2) * Cue_Pred.coeff(count,2)) + (dataset(iTrial,9) * Cue_Pred.coeff(count,3))));
            end
            count = count + 1;
        end
    end
end

%%
figure
for iPlot = 1:10
    subplot(2,5,iPlot)
    scatter(1:length(Cue_Pred.pred{iPlot}),Cue_Pred.pred{iPlot}(:,1))
    hold on;
    scatter(1:length(Cue_Pred.pred{iPlot}),Cue_Pred.pred{iPlot}(:,2),'r')
    scatter(1:length(Cue_Pred.pred{iPlot}),Cue_Pred.pred{iPlot}(:,3),'g')
    ylim([-.5 1.5]); %title('Trial Length');
    % set(gca,'TickLength', [0 0]); box off;
    set(gca,'XTickLabel',[])
    %     ylabel('Trial length (s)');  xlabel('Cue type'); %'XTickLabelRotation',90,
    box off;
    h = gca;
    h.XRuler.TickLength = 0;
    set(h,'FontSize',18);
end

figure
for iPlot = 1:2
    subplot(2,4,iPlot)
    scatter(1:length(Cue_Pred.pred{iPlot+8}),Cue_Pred.pred{iPlot+8}(:,1))
    hold on;
    scatter(1:length(Cue_Pred.pred{iPlot+8}),Cue_Pred.pred{iPlot+8}(:,2),'r')
    scatter(1:length(Cue_Pred.pred{iPlot+8}),Cue_Pred.pred{iPlot+8}(:,3),'g')
    ylim([-.5 1.5]); %title('Trial Length');
    % set(gca,'TickLength', [0 0]); box off;
    set(gca,'XTickLabel',[])
    %     ylabel('Trial length (s)');  xlabel('Cue type'); %'XTickLabelRotation',90,
    box off;
    h = gca;
    h.XRuler.TickLength = 0;
    set(h,'FontSize',18);
end
%%
count_crit = 1;
for iCell = 1:length(mdl)
    MdlComp{t_count}.comp(iCell,1) = mdl{iCell}.ALL.OneBack.ModelCriterion.AIC - mdl{iCell}.ALL.OneFRate.ModelCriterion.AIC;
    MdlComp{t_count}.comp(iCell,2) =  mdl{iCell}.ALL.OneBack.ModelCriterion.BIC - mdl{iCell}.ALL.OneFRate.ModelCriterion.BIC;
    MdlComp{t_count}.comp(iCell,3) = mdl{iCell}.ALL.OneBack.Deviance - mdl{iCell}.ALL.OneFRate.Deviance;
    MdlComp{t_count}.comp(iCell,4) = mdl{iCell}.ALL.OneBack.Rsquared.Adjusted - mdl{iCell}.ALL.OneFRate.Rsquared.Adjusted;
    if MdlComp{t_count}.ALL.FRate(iCell,5) == 1
        sigComp(count_crit,:) = MdlComp{t_count}.comp(iCell,:);
        count_crit = count_crit + 1;
    end
end

%%
figure
hold on;
for iMdl = 1:3
    scatter(ones(1,443)*iMdl,MdlComp{t_count}.comp(:,iMdl),'r')
    plot([iMdl-.30 iMdl+.30],[mean(MdlComp{t_count}.comp(:,iMdl)) mean(MdlComp{t_count}.comp(:,iMdl))],'k')
    scatter(ones(1,length(sigComp))*iMdl,sigComp(:,iMdl),'c')
    plot([iMdl-.30 iMdl+.30],[mean(sigComp(:,iMdl)) mean(sigComp(:,iMdl))],'b')
end

plot(0:.05:4,0,'.','color','black')

xlim([0 4]); ylim([-20 20]); title('Improvement to model fit from addition of FRate to 1-back');
set(gca,'XTickLabel',{'','AIC','BIC','Deviance'});
ylabel('Improvement to model fit (AU)'); xlabel('Model criterion');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%%
figure
hold on;

scatter(ones(1,443)*1,MdlComp{t_count}.ALL.predicted(:,3),'r')
plot([1-.30 1+.30],[mean(MdlComp{t_count}.ALL.predicted(:,3)) mean(MdlComp{t_count}.ALL.predicted(:,3))],'k')
scatter(ones(1,12)*2,FRate_cells2,'r')
plot([2-.30 2+.30],[mean(FRate_cells2) mean(FRate_cells2)],'k')
scatter(ones(1,443)*3,MdlComp{t_count}.ALL.predicted(:,2),'r')
plot([3-.30 3+.30],[mean(MdlComp{t_count}.ALL.predicted(:,2)) mean(MdlComp{t_count}.ALL.predicted(:,2))],'k')
scatter(ones(1,15)*4,FRate_cells,'r')
plot([4-.30 4+.30],[mean(FRate_cells) mean(FRate_cells)],'k')

plot(0:.05:10,.5,'.','color','black')

xlim([0 5]); ylim([0 1]); title('Proportion of correctly predicted next trial reward type');
set(gca,'XTickLabel',{'','1-back' 'FRate SIG' 'FRate all','FRate SIG'});
ylabel('Proportion correctly predicted'); xlabel('Model specs');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%%


%%
figure
hold on;
for iMdl = 1:length(mdl_identifier)
    scatter(ones(1,443)*iMdl,MdlComp{t_count}.ALL.predicted(:,iMdl),'r')
    plot([iMdl-.30 iMdl+.30],[mean(MdlComp{t_count}.ALL.predicted(:,iMdl)) mean(MdlComp{t_count}.ALL.predicted(:,iMdl))],'k')
    scatter(ones(1,length(sigF.ALL))*iMdl,sigF.ALL(:,iMdl),'c')
    plot([iMdl-.30 iMdl+.30],[mean(sigF.ALL(:,iMdl)) mean(sigF.ALL(:,iMdl))],'b')
end

plot(0:.05:10,.5,'.','color','black')

xlim([0 9]); ylim([0 1]); title('Proportion of correctly predicted next trial reward type');
set(gca,'XTickLabel',{'','1-back','2-back','App','RRate','FRate','1-FRate','App-FRate','R-FRate','Int'});
ylabel('Proportion correctly predicted'); xlabel('Model specs');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%%
%%
subplot(1,2,2)
hold on;
for iMdl = 1:length(mdl_identifier)
    scatter(ones(1,length(sigF.ALL))*iMdl,sigF.ALL(:,iMdl),'r')
    plot([iMdl-.30 iMdl+.30],[mean(sigF.ALL(:,iMdl)) mean(sigF.ALL(:,iMdl))],'k')
end

plot(0:.05:10,.5,'.','color','black')

xlim([0 9]); ylim([0 1]); title('Proportion of correctly predicted next trial reward type');
set(gca,'XTickLabel',{'','1-back','2-back','App','RRate','FRate','1-FRate','App-FRate','R-FRate','Int'});
ylabel('Proportion correctly predicted'); xlabel('Model specs');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%%
for iTime = 1:11
    FRate_count(iTime) = MdlComp{iTime}.ALL.Count;
    FRate_count_cue(iTime) = MdlComp{iTime}.cuemod.Count;
end

%%
figure
subplot(1,2,1)
plot(-500:100:500,FRate_count/433)
xlim([-600 600]); ylim([0 .5]); title('FRate sig predictor (all 443 units)');
xlabel('Start of 500 ms bin');
ylabel('Proportion of units');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

subplot(1,2,2)
plot(-500:100:500,FRate_count_cue/133)
xlim([-600 600]); ylim([0 .5]); title('FRate sig predicotr (cue modulated 133 units)');
xlabel('Start of 500 ms bin');
ylabel('Proportion of units');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%% Filtering sig predicts through predictability
mat_files = dir('*.mat');
t_count = 0;
modelspec = {'Outcome ~ 1 + Previous','Outcome ~ 1 + Previous + FiringRate','Outcome ~ 1 + FiringRate'};
mdl_identifier = {'OneBack','OneFRate','FRate'};

for iTime = -.5:.1:.5
    load(cat(2,'E:\Jimmie\Jimmie\Analysis\2018-06-06-GLM_cueon_',num2str(iTime),'-replication.mat'));
    
    t_count = t_count + 1;
    
    for iMdl = 1:length(mdl_identifier)
        count{iMdl} = 1;
        count_bdrift{iMdl} = 1;
        count_cuemod{iMdl} = 1;
    end
    
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    for iCell = 1:length(ALL_matrix)
        if ALL_matrix(iCell,3) == 1
            
            load(mat_files(iCell).name);
            mat_overview.fname{iCell} = mat_files(iCell).name;
            disp(cat(2,num2str(iCell),'/',num2str(length(dir('*.mat')))));
            new_v_old = strcmp(mat_overview.fname{iCell}(1:4),'R060');
            
            clear dataset dataset2 ds %mdl
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo_block1.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                            
                        else
                            dataset(trials_count,1) = metadata.TrialInfo_block1.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo_block1.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                        switch sesh.block_order %modality
                            case 1
                                dataset(trials_count,4) = 1;
                            case 2
                                dataset(trials_count,4) = 2;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo_block1.summary(jj,5); %arm
                        dataset(trials_count,6) = metadata.TrialInfo_block1.approached(jj); %behav - app
                        dataset(trials_count,7) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                        dataset(trials_count,8) = jj;
                        dataset(trials_count,9) = firing_rate(trials_count); %FR
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo_block2.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo_block2.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo_block2.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo_block2.rewarded(jj);
                        switch sesh.block_order
                            case 1
                                dataset(trials_count,4) = 2;
                            case 2
                                dataset(trials_count,4) = 1;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo_block2.summary(jj,5);
                        dataset(trials_count,6) = metadata.TrialInfo_block2.approached(jj);
                        dataset(trials_count,7) = metadata.TrialInfo_block2.summary(jj,15);
                        dataset(trials_count,8) = b1_length+jj;
                        dataset(trials_count,9) = firing_rate(trials_count);
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo{1,1}.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo{1,1}.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo{1,1}.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                        switch sesh.block_order %modality
                            case 1
                                dataset(trials_count,4) = 1;
                            case 2
                                dataset(trials_count,4) = 2;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                        dataset(trials_count,6) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                        dataset(trials_count,7) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                        dataset(trials_count,8) = jj;
                        dataset(trials_count,9) = firing_rate(trials_count); %FR
                        
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
                            dataset(trials_count,2) = NaN; %prev trial
                            dataset(trials_count,10) = NaN;
                            dataset(trials_count,11) = NaN;
                        elseif jj == 2
                            dataset(trials_count,1) = NaN; %prev trial
                            dataset(trials_count,2) =  metadata.TrialInfo{1,2}.rewarded(jj-1); %prev trial
                            dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2);
                        else
                            dataset(trials_count,1) = metadata.TrialInfo{1,2}.rewarded(jj-2);
                            dataset(trials_count,2) = metadata.TrialInfo{1,2}.rewarded(jj-1);
                            dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                            dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                        end
                        
                        dataset(trials_count,3) = metadata.TrialInfo{1,2}.rewarded(jj);
                        switch sesh.block_order
                            case 1
                                dataset(trials_count,4) = 2;
                            case 2
                                dataset(trials_count,4) = 1;
                        end
                        dataset(trials_count,5) = metadata.TrialInfo{1,2}.summary(jj,5);
                        dataset(trials_count,6) = metadata.TrialInfo{1,2}.approached(jj);
                        dataset(trials_count,7) = metadata.TrialInfo{1,2}.summary(jj,15);
                        dataset(trials_count,8) = b1_length+jj;
                        dataset(trials_count,9) = firing_rate(trials_count);
                        
                        trials_count = trials_count + 1;
                        
                    end
                    
            end
            
            %         ds = mat2dataset(dataset,'VarNames',{'PrevTwo','Previous','Outcome','Modality','Location','Approach','Latency','Trial','FiringRate','PrevApp','RRate'});
            dataset2 = cat(2,dataset(:,2),dataset(:,9),dataset(:,3));
            ds = mat2dataset(dataset2,'VarNames',{'Previous','FiringRate','Outcome'});
            
            %         mdl_comp = mnrfit
            for iMdl = 1:length(mdl_identifier)
                
                mdl_comp{t_count,iCell}.ALL.(mdl_identifier{iMdl}) = fitglm(ds,modelspec{iMdl},'Link','logit','Distribution','binomial');
                %             MdlComp{t_count}.ALL.Deviance(iCell,iMdl) =  mdl{t_count,iCell}.ALL.(mdl_identifier{iMdl}).Deviance;
                %             MdlComp{t_count}.ALL.pvalue(iCell,iMdl) = 1-chi2cdf(MdlComp{t_count}.ALL.Deviance(iCell,1)-MdlComp{t_count}.ALL.Deviance(iCell,iMdl),1);
                out_og = mdl_comp{t_count,iCell}.ALL.(mdl_identifier{iMdl}).Fitted;
                out = categorical(out_og.Response > 0.5);
                data = categorical(dataset(:,3) > 0.5);
                residuals{iMdl} = abs(dataset(:,3) - out_og.Response);
                residual = nanmean(abs(dataset(:,3) - out_og.Response));
                residual_sq = nansum((dataset(:,3) - out_og.Response).^2);
                MdlComp{t_count}.ALL.predicted(iCell,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
                MdlComp{t_count}.ALL.predRatio(iCell,iMdl) = MdlComp{t_count}.ALL.predicted(iCell,iMdl) / MdlComp{t_count}.ALL.predicted(iCell,1);
                MdlComp{t_count}.ALL.residual(iCell,iMdl) = residual;
                MdlComp{t_count}.ALL.residualSS(iCell,iMdl) = residual_sq;
                MdlComp{t_count}.ALL.residualsNull(iCell,iMdl) = nanmean(residuals{iMdl} - residuals{1});
                %             for iCoeff = 1:length(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames)
                %                 if strcmp(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FRate') == 1
                %                     if mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                %                         MdlComp{t_count}.ALL.FRate(iSesh,iMdl) = 1;
                %                         MdlComp{t_count}.ALL.Count = count{iMdl};
                %                         MdlComp{t_count}.ALL.FRate_cells{iMdl}(count{iMdl}) = MdlComp{t_count}.ALL.predicted(iSesh,iMdl);
                %                         count{iMdl} = count{iMdl} + 1;
                %                     end
                %                 end
                %             end
                
            end
        end
    end
end

%%
t_count = 0;
for iTime = -.5:.1:.5
    t_count = t_count + 1;
    diff_count = 1;
    for iCell = 1:length(MdlComp{t_count}.ALL.predicted)
        if MdlComp{t_count}.ALL.predicted(iCell,1) > 0
            MdlDiff{t_count}(diff_count) = MdlComp{t_count}.ALL.predicted(iCell,2) - MdlComp{t_count}.ALL.predicted(iCell,1);
            
            MdlFRate{t_count}(diff_count) = MdlComp{t_count}.ALL.predicted(iCell,3);
            diff_count = diff_count + 1;
        end
    end
    MdlAvg(t_count) = mean(MdlDiff{t_count});
    MdlSEM(t_count) = std(MdlDiff{t_count})/  sqrt(numel(MdlDiff{t_count}));
    MdlAvgFR(t_count) = mean(MdlFRate{t_count});
    MdlSEMFR(t_count) = std(MdlFRate{t_count})/  sqrt(numel(MdlFRate{t_count}));
end

%%
figure
hold on;
% for t_count = 1:11
%     scatter(ones(1,length(MdlDiff{t_count}))*t_count,MdlDiff{t_count},'r')
%     plot([t_count-.30 t_count+.30],[mean(MdlDiff{t_count}) mean(MdlDiff{t_count})],'k')
% end

for t_count = 1:11
    scatter(ones(1,length(MdlFRate{t_count}))*t_count,MdlFRate{t_count},'r')
    plot([t_count-.30 t_count+.30],[mean(MdlFRate{t_count}) mean(MdlFRate{t_count})],'k')
end

shadedErrorBar(1:11,MdlAvgFR,MdlSEMFR,'-b',1);

plot(0:.05:11,.5,'.','color','black')

xlim([0 11]); ylim([0.2 .8]); title('Proportion of correctly predicted next trial reward type');
set(gca,'XTickLabel',{'','1-back','2-back','App','RRate','FRate','1-FRate','App-FRate','R-FRate','Int'});
ylabel('Proportion correctly predicted'); xlabel('Model specs');
box off;
h = gca;
h.XRuler.TickLength = 0;
set(h,'FontSize',18);

%%

figure
shadedErrorBar(-.25:.1:.75,MdlAvgFR,MdlSEMFR,'-b',1);

%% cross val
rng('shuffle')
warning('error', 'stats:glmfit:IterationLimit'); % turn GLM warning into error so that you can use try catch to skip to next iteration of loop.

mat_files = dir('*.mat');
t_count = 0;
mdl_identifier = {'OneBack','OneFRate'};

for iTime = -.5:.1:.5
    t_count = t_count + 1;
    count = 1;
    
    time_window_start = iTime; %starting time window for analysis, 0 = time zero
    time_window_end = time_window_start + .5;
    epoch_start = 0;
    
    for iCell = 1:length(mat_files)
        
        load(mat_files(iCell).name);
        mat_overview.fname{iCell} = mat_files(iCell).name;
        disp(cat(2,num2str(iCell),'/',num2str(length(dir('*.mat'))),' (',num2str(iTime),' epoch)'));
        new_v_old = strcmp(mat_overview.fname{iCell}(1:4),'R060');
        
        clear dataset x y %mdl
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo_block1.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                        
                    else
                        dataset(trials_count,1) = metadata.TrialInfo_block1.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo_block1.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo_block1.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo_block1.rewarded(jj); %outcome
                    switch sesh.block_order %modality
                        case 1
                            dataset(trials_count,4) = 1;
                        case 2
                            dataset(trials_count,4) = 2;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo_block1.summary(jj,5); %arm
                    dataset(trials_count,6) = metadata.TrialInfo_block1.approached(jj); %behav - app
                    dataset(trials_count,7) = metadata.TrialInfo_block1.summary(jj,15); %behav - trial length
                    dataset(trials_count,8) = jj;
                    dataset(trials_count,9) = firing_rate(trials_count); %FR
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo_block2.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo_block2.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo_block2.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo_block2.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo_block2.rewarded(jj);
                    switch sesh.block_order
                        case 1
                            dataset(trials_count,4) = 2;
                        case 2
                            dataset(trials_count,4) = 1;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo_block2.summary(jj,5);
                    dataset(trials_count,6) = metadata.TrialInfo_block2.approached(jj);
                    dataset(trials_count,7) = metadata.TrialInfo_block2.summary(jj,15);
                    dataset(trials_count,8) = b1_length+jj;
                    dataset(trials_count,9) = firing_rate(trials_count);
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo{1,1}.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo{1,1}.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo{1,1}.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo{1,1}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo{1,1}.rewarded(jj); %outcome
                    switch sesh.block_order %modality
                        case 1
                            dataset(trials_count,4) = 1;
                        case 2
                            dataset(trials_count,4) = 2;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo{1,1}.summary(jj,5); %arm
                    dataset(trials_count,6) = metadata.TrialInfo{1,1}.approached(jj); %behav - app
                    dataset(trials_count,7) = metadata.TrialInfo{1,1}.summary(jj,15); %behav - trial length
                    dataset(trials_count,8) = jj;
                    dataset(trials_count,9) = firing_rate(trials_count); %FR
                    
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
                        dataset(trials_count,2) = NaN; %prev trial
                        dataset(trials_count,10) = NaN;
                        dataset(trials_count,11) = NaN;
                    elseif jj == 2
                        dataset(trials_count,1) = NaN; %prev trial
                        dataset(trials_count,2) =  metadata.TrialInfo{1,2}.rewarded(jj-1); %prev trial
                        dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2);
                    else
                        dataset(trials_count,1) = metadata.TrialInfo{1,2}.rewarded(jj-2);
                        dataset(trials_count,2) = metadata.TrialInfo{1,2}.rewarded(jj-1);
                        dataset(trials_count,10) =  metadata.TrialInfo{1,2}.approached(jj-1);
                        dataset(trials_count,11) = dataset(trials_count,2) + dataset(trials_count,1);
                    end
                    
                    dataset(trials_count,3) = metadata.TrialInfo{1,2}.rewarded(jj);
                    switch sesh.block_order
                        case 1
                            dataset(trials_count,4) = 2;
                        case 2
                            dataset(trials_count,4) = 1;
                    end
                    dataset(trials_count,5) = metadata.TrialInfo{1,2}.summary(jj,5);
                    dataset(trials_count,6) = metadata.TrialInfo{1,2}.approached(jj);
                    dataset(trials_count,7) = metadata.TrialInfo{1,2}.summary(jj,15);
                    dataset(trials_count,8) = b1_length+jj;
                    dataset(trials_count,9) = firing_rate(trials_count);
                    
                    trials_count = trials_count + 1;
                    
                end
                
        end
        
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
                if RANK.two.Trial > 975 || RANK.two.Trial < 26
                    if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                        
%                         shuff_FRate = [];
%                         for iShuff = 1:100
%                             shuff_FRate(:,1) = datasample(dataset(:,9),length(dataset(:,9)),'Replace',false);
                            x{1} = dataset(:,2);
                            x{2} = cat(2,dataset(:,2),dataset(:,9));
%                             x{2} = cat(2,dataset(:,2),shuff_FRate(:,1));
                            y = dataset(:,3);
                            
                            for iMdl = 1:length(mdl_identifier)
                                try
                                    yfit = @(xtrain,ytrain,xtest)(glmCUE(xtrain,ytrain,xtest)');
                                    cvMcr{t_count}(count,iMdl) = crossval('mcr',x{iMdl},y,'predfun',yfit);
%                                cvMse.(cat(2,'s',num2str(iShuff))){t_count}(count,iMdl) = crossval('mse',x{iMdl},y,'predfun',yfit);
                                catch
                                    cvMcr{t_count}(count,iMdl) = NaN;
%                                      cvMse.(cat(2,'s',num2str(iShuff))){t_count}(count,iMdl) = NaN;
                                end
                            end
%                         end
                        count = count + 1;
                    end
                end
        end
    end
%     for iShuff = 1:100
        cvMcr_diff(:,t_count) = cvMcr{t_count}(:,1) - cvMcr{t_count}(:,2);
%         cvMse_diff.(cat(2,'s',num2str(iShuff)))(:,t_count) = cvMse.(cat(2,'s',num2str(iShuff))){t_count}(:,1) - cvMse.(cat(2,'s',num2str(iShuff))){t_count}(:,2);
%     end
end
%%
for t_count = 1:11
 cvMse_diff(:,t_count) = cvMse{t_count}(:,1) - cvMse{t_count}(:,2);
end

%%
for iShuff = 1:100
    cvMse_SHUFF.ALL(iShuff,:) = nanmean(cvMse_diff_SHUFF.(cat(2,'s',num2str(iShuff))));
end

cvMse_SHUFF.MEAN = mean(cvMse_SHUFF.ALL);
cvMse_SHUFF.SEM = std(cvMse_SHUFF.ALL) / sqrt(length(cvMse_SHUFF.ALL));

%%

figure
violin(cvMse_SHUFF.ALL,'FaceColor',[0 0 0],'LineColor',[1 0 0])
hold on;
violin(cvMse_diff);
set(gca,'XTickLabel',{'Pre-cue',' ',' ',' ',' ','Cue',' ',' ',' ',' ',' '})
% ylim([0 12]);
ylabel('Decrease in MSE')
title('Comparison of 1-back and FRate with FRate')
set(gca,'FontSize',18);