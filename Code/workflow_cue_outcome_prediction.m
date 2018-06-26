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
            mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9});
            MdlComp.ALL.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Deviance;
            MdlComp.ALL.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp.ALL.Deviance(iSesh,1)-MdlComp.ALL.Deviance(iSesh,iMdl),1);
            out = mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Fitted;
            out = categorical(out.Response > 0.5);
            data = categorical(dataset(:,3) > 0.5);
            MdlComp.ALL.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
            MdlComp.ALL.predRatio(iSesh,iMdl) = MdlComp.ALL.predicted(iSesh,iMdl) / MdlComp.ALL.predicted(iSesh,1);
            for iCoeff = 1:length(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames)
                if strcmp(mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FiringRate') == 1
                    if mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                        MdlComp.ALL.FRate(iSesh,iMdl) = 1;
                        MdlComp.ALL.FRate_cells{iMdl}(count{iMdl}) = MdlComp.ALL.predicted(iSesh,iMdl);
                        count{iMdl} = count{iMdl} + 1;
                    else
                        MdlComp.ALL.FRate(iSesh,iMdl) = 0;
                    end
                end
            end
            
            if bdrift == 1
                mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9});
                MdlComp.bdrift.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Deviance;
                MdlComp.bdrift.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp.bdrift.Deviance(iSesh,1)-MdlComp.bdrift.Deviance(iSesh,iMdl),1);
                out = mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Fitted;
                out = categorical(out.Response > 0.5);
                data = categorical(dataset(:,3) > 0.5);
                MdlComp.bdrift.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
                MdlComp.bdrift.predRatio(iSesh,iMdl) = MdlComp.bdrift.predicted(iSesh,iMdl) / MdlComp.bdrift.predicted(iSesh,1);
                for iCoeff = 1:length(mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).CoefficientNames)
                    if strcmp(mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FiringRate') == 1
                        if mdl{t_count,iSesh}.bdrift.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                            MdlComp.bdrift.FRate(iSesh,iMdl) = 1;
                                                    MdlComp.bdrift.FRate_cells{iMdl}(count_bdrift{iMdl}) = MdlComp.bdrift.predicted(iSesh,iMdl);
                        count_bdrift{iMdl} = count_bdrift{iMdl} + 1;
                        else
                            MdlComp.bdrift.FRate(iSesh,iMdl) = 0;
                        end
                    end
                end
                
                if cuemod == 1
                mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}) = stepwiseglm(ds,modelspec{9});
                MdlComp.cuemod.Deviance(iSesh,iMdl) =  mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Deviance;
                MdlComp.cuemod.pvalue(iSesh,iMdl) = 1-chi2cdf(MdlComp.cuemod.Deviance(iSesh,1)-MdlComp.cuemod.Deviance(iSesh,iMdl),1);
                out = mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Fitted;
                out = categorical(out.Response > 0.5);
                data = categorical(dataset(:,3) > 0.5);
                MdlComp.cuemod.predicted(iSesh,iMdl) = sum(out == data)/length(dataset(:,3)); % count number of correct predictions
                MdlComp.cuemod.predRatio(iSesh,iMdl) = MdlComp.cuemod.predicted(iSesh,iMdl) / MdlComp.cuemod.predicted(iSesh,1);
                for iCoeff = 1:length(mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).CoefficientNames)
                    if strcmp(mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).CoefficientNames(iCoeff),'FiringRate') == 1
                        if mdl{t_count,iSesh}.cuemod.(mdl_identifier{iMdl}).Coefficients{iCoeff,4} < .01
                            MdlComp.cuemod.FRate(iSesh,iMdl) = 1;
                                                    MdlComp.cuemod.FRate_cells{iMdl}(count_cuemod{iMdl}) = MdlComp.cuemod.predicted(iSesh,iMdl);
                        count_cuemod{iMdl} = count_cuemod{iMdl} + 1;
                        else
                            MdlComp.cuemod.FRate(iSesh,iMdl) = 0;
                        end
                    end
                end
            end
            end
            
            
                        catch
                            mdl{t_count,iSesh}.ALL.(mdl_identifier{iMdl}) = [];
                        end
%         end
        
        %         MdlComp.ALL.Proportions(iSesh,1) = sum(dataset(:,2) == 0 & dataset(:,3) == 0);
        %         MdlComp.ALL.Proportions(iSesh,2) = sum(dataset(:,2) == 1 & dataset(:,3) == 0);
        %         MdlComp.ALL.Proportions(iSesh,3) = sum(dataset(:,2) == 0 & dataset(:,3) == 1);
        %         MdlComp.ALL.Proportions(iSesh,4) = sum(dataset(:,2) == 1 & dataset(:,3) == 1);
        %         MdlComp.ALL.CondP(iSesh,1) = MdlComp.ALL.Proportions(iSesh,1) / (MdlComp.ALL.Proportions(iSesh,1) + MdlComp.ALL.Proportions(iSesh,2));
        %         MdlComp.ALL.CondP(iSesh,2) = MdlComp.ALL.Proportions(iSesh,2) / (MdlComp.ALL.Proportions(iSesh,1) + MdlComp.ALL.Proportions(iSesh,2));
        %         MdlComp.ALL.CondP(iSesh,3) = MdlComp.ALL.Proportions(iSesh,3) / (MdlComp.ALL.Proportions(iSesh,3) + MdlComp.ALL.Proportions(iSesh,4));
        %         MdlComp.ALL.CondP(iSesh,4) = MdlComp.ALL.Proportions(iSesh,4) / (MdlComp.ALL.Proportions(iSesh,3) + MdlComp.ALL.Proportions(iSesh,4));
    end
end

% ds.Outcome = categorical(ds.Outcome);
% ds.Previous = categorical(ds.Previous);
% modelspec = 'Outcome ~ 1 + Previous';
% glm = fitglm(ds,modelspec);
%
% p = 1-chi2cdf(glm_linear.Deviance-glm_linear2.Deviance,1)

%%
MdlComp.ALL.FRate_count = 0;

for iSesh = 1:length(mdl)
    for iCoeff = 1:length(mdl{t_count,iSesh}.ALL.OneBack.CoefficientNames)
                    if strcmp(mdl{t_count,iSesh}.ALL.OneBack.CoefficientNames(iCoeff),'FRate') == 1
                        if mdl{t_count,iSesh}.ALL.OneBack.Coefficients{iCoeff,4} < .01
                            MdlComp.ALL.FRate(iSesh) = 1;
                                                    MdlComp.ALL.FRate_count = MdlComp.ALL.FRate_count + 1;
                        else
                            MdlComp.ALL.FRate(iSesh) = 0;
                        end
                    end
                end
end

%%
MdlComp.bdrift.FRate_count = 0;

for iSesh = 1:length(mdl)
    if isfield(mdl{t_count,iSesh},'bdrift') == 1
    for iCoeff = 1:length(mdl{t_count,iSesh}.bdrift.OneBack.CoefficientNames)
                    if strcmp(mdl{t_count,iSesh}.bdrift.OneBack.CoefficientNames(iCoeff),'FRate') == 1
                        if mdl{t_count,iSesh}.bdrift.OneBack.Coefficients{iCoeff,4} < .01
                            MdlComp.bdrift.FRate(iSesh) = 1;
                                                    MdlComp.bdrift.FRate_count = MdlComp.bdrift.FRate_count + 1;
                        else
                            MdlComp.bdrift.FRate(iSesh) = 0;
                        end
                    end
    end
    end
end

%%
for iMdl = 1:length(MdlComp.ALL.FRate(1,:))
MdlComp.ALL.SummaryRAW.FRate(iMdl) = sum(MdlComp.ALL.FRate(:,iMdl) == 1);
MdlComp.bdrift.SummaryRAW.FRate(iMdl) = sum(MdlComp.bdrift.FRate(:,iMdl) == 1);
MdlComp.cuemod.SummaryRAW.FRate(iMdl) = sum(MdlComp.cuemod.FRate(:,iMdl) == 1);

for iMdl2 = 1:length(MdlComp.ALL.FRate(1,:))
    MdlComp.ALL.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp.ALL.FRate(:,iMdl) == 1 & MdlComp.ALL.FRate(:,iMdl2) == 1);
        MdlComp.ALL.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp.ALL.FRate(:,iMdl) == 1 & MdlComp.ALL.FRate(:,iMdl2) == 0);
        
        MdlComp.bdrift.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp.bdrift.FRate(:,iMdl) == 1 & MdlComp.bdrift.FRate(:,iMdl2) == 1);
        MdlComp.bdrift.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp.bdrift.FRate(:,iMdl) == 1 & MdlComp.bdrift.FRate(:,iMdl2) == 0);
        
        MdlComp.cuemod.SummaryRAW.Combined(iMdl,iMdl2) = sum(MdlComp.cuemod.FRate(:,iMdl) == 1 & MdlComp.cuemod.FRate(:,iMdl2) == 1);
        MdlComp.cuemod.SummaryRAW.Exclusive(iMdl,iMdl2) = sum(MdlComp.cuemod.FRate(:,iMdl) == 1 & MdlComp.cuemod.FRate(:,iMdl2) == 0);
end
end



%%
count_all = 1;
count_bd = 1;
count_cm = 1;

for iCell = 1:442%length(mdl)
    if MdlComp.ALL.FRate(iCell,5) == 1
        sigF.ALL(count_all,:) = MdlComp.ALL.predicted(iCell,:);
        count_all = count_all + 1;
    end
    if MdlComp.bdrift.FRate(iCell,5) == 1
        sigF.bd(count_bd,:) = MdlComp.bdrift.predicted(iCell,:);
        count_bd = count_bd + 1;
    end
    if MdlComp.cuemod.FRate(iCell,5) == 1
        sigF.cm(count_cm,:) = MdlComp.cuemod.predicted(iCell,:);
        count_cm = count_cm + 1;
    end
end
    
% MdlComp.ALL.SummaryRAW.Combined = sum(MdlComp.ALL.FRate(:,3) == 1);
% MdlComp.ALL.SummaryRAW.FRateOnly = sum(MdlComp.ALL.FRate(:,2) == 1 & MdlComp.ALL.FRate(:,3) == 0);
% MdlComp.ALL.SummaryRAW.CombinedOnly = sum(MdlComp.ALL.FRate(:,3) == 1 & MdlComp.ALL.FRate(:,2) == 0);
% MdlComp.ALL.SummaryRAW.Both = sum(MdlComp.ALL.FRate(:,2) == 1 & MdlComp.ALL.FRate(:,3) == 1);
% MdlComp.ALL.SummaryPROP.FRate = sum(MdlComp.ALL.FRate(:,2) == 1) /133;
% MdlComp.ALL.SummaryPROP.Combined = sum(MdlComp.ALL.FRate(:,3) == 1) /133;
% MdlComp.ALL.SummaryPROP.FRateOnly = sum(MdlComp.ALL.FRate(:,2) == 1 & MdlComp.ALL.FRate(:,3) == 0) /133;
% MdlComp.ALL.SummaryPROP.CombinedOnly = sum(MdlComp.ALL.FRate(:,3) == 1 & MdlComp.ALL.FRate(:,2) == 0) /133;
% MdlComp.ALL.SummaryPROP.Both = sum(MdlComp.ALL.FRate(:,2) == 1 & MdlComp.ALL.FRate(:,3) == 1) /133;
%%

% count = 1;
% for iCell = 1:length(MdlComp.ALL.FRate)
%     if MdlComp.ALL.FRate(iCell,2) == 1
%         FRate_cells(count) = MdlComp.ALL.predicted(iCell,2);
%         count = count + 1;
%     end
% end
% 
% count = 1;
% for iCell = 1:length(MdlComp.ALL.FRate)
%     if MdlComp.ALL.FRate(iCell,3) == 1
%         FRate_cells2(count) = MdlComp.ALL.predicted(iCell,3);
%         count = count + 1;
%     end
% end
%%
count_crit = 1;
for iCell = 1:length(mdl)
    MdlComp.comp(iCell,1) = mdl{iCell}.ALL.OneBack.ModelCriterion.AIC - mdl{iCell}.ALL.OneFRate.ModelCriterion.AIC;
    MdlComp.comp(iCell,2) =  mdl{iCell}.ALL.OneBack.ModelCriterion.BIC - mdl{iCell}.ALL.OneFRate.ModelCriterion.BIC;
    MdlComp.comp(iCell,3) = mdl{iCell}.ALL.OneBack.Deviance - mdl{iCell}.ALL.OneFRate.Deviance;
    MdlComp.comp(iCell,4) = mdl{iCell}.ALL.OneBack.Rsquared.Adjusted - mdl{iCell}.ALL.OneFRate.Rsquared.Adjusted;
    if MdlComp.ALL.FRate(iCell,5) == 1
        sigComp(count_crit,:) = MdlComp.comp(iCell,:);
        count_crit = count_crit + 1;
    end
end

%%
figure
hold on;
for iMdl = 1:3
    scatter(ones(1,443)*iMdl,MdlComp.comp(:,iMdl),'r')
    plot([iMdl-.30 iMdl+.30],[mean(MdlComp.comp(:,iMdl)) mean(MdlComp.comp(:,iMdl))],'k')
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

scatter(ones(1,443)*1,MdlComp.ALL.predicted(:,3),'r')
plot([1-.30 1+.30],[mean(MdlComp.ALL.predicted(:,3)) mean(MdlComp.ALL.predicted(:,3))],'k')
scatter(ones(1,12)*2,FRate_cells2,'r')
plot([2-.30 2+.30],[mean(FRate_cells2) mean(FRate_cells2)],'k')
scatter(ones(1,443)*3,MdlComp.ALL.predicted(:,2),'r')
plot([3-.30 3+.30],[mean(MdlComp.ALL.predicted(:,2)) mean(MdlComp.ALL.predicted(:,2))],'k')
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
    scatter(ones(1,443)*iMdl,MdlComp.ALL.predicted(:,iMdl),'r')
    plot([iMdl-.30 iMdl+.30],[mean(MdlComp.ALL.predicted(:,iMdl)) mean(MdlComp.ALL.predicted(:,iMdl))],'k')
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