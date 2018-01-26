count = 1;
cells_to_use = [];
for iCell = 1:length(block_drift.summary)
    if block_drift.summary(iCell) == 0
        cells_to_use(count) = iCell;
        count = count + 1;
    end
end

for iShuff = 1:1000
disp(iShuff)
%% average for modality
num_cells = length(find(ALL_matrix(:,1) == 1));
rng('shuffle')
load_vars = datasample(cells_to_use,num_cells,'Replace',false);
% load_vars = randperm(length(ALL_matrix),num_cells);
% load_vars = find(ALL_matrix(:,1) == 1);

rule_count = 1;
rule_count_INC = 1;
rule_count_DEC = 1;

mat_files = dir('*.mat');
for kk = 1:length(load_vars)
    load(mat_files(load_vars(kk)).name);
    mat_overview.fname{kk} = mat_files(load_vars(kk)).name;
%     disp(cat(2,num2str(kk),'/',num2str(length(load_vars))));   

                temp_block1_MEAN = [];
                temp_block1_SEM = [];
                temp_block2_MEAN = [];
                temp_block2_SEM = [];
                
                rule_encoding.zscore.MEAN(rule_count) = mean(PETH.Trial.MEAN.trials_PETH);
                rule_encoding.zscore.STD(rule_count) = std(PETH.Trial.MEAN.trials_PETH);
                               
                switch sesh.block_order
                    case 1
                        temp_block1_MEAN = PETH.Trial.MEAN.trials_light_PETH;
                        temp_block1_SEM = PETH.Trial.SEM.trials_light_PETH;
                        temp_block2_MEAN = PETH.Trial.MEAN.trials_sound_PETH;
                        temp_block2_SEM = PETH.Trial.SEM.trials_sound_PETH;
                    case 2
                        temp_block1_MEAN = PETH.Trial.MEAN.trials_sound_PETH;
                        temp_block1_SEM = PETH.Trial.SEM.trials_sound_PETH;
                        temp_block2_MEAN = PETH.Trial.MEAN.trials_light_PETH;
                        temp_block2_SEM = PETH.Trial.SEM.trials_light_PETH;
                end
                
                if mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                    if mean(FRATE.Cue.Trial_firing_rate_block1) > mean(FRATE.Cue.Trial_firing_rate_block2)
                        rule_encoding.RAW.pref.ALL.MEAN(:,rule_count) = temp_block1_MEAN;
                        rule_encoding.RAW.pref.ALL.SEM(:,rule_count) = temp_block1_SEM;
                        rule_encoding.RAW.nonpref.ALL.MEAN(:,rule_count) = temp_block2_MEAN;
                        rule_encoding.RAW.nonpref.ALL.SEM(:,rule_count) = temp_block2_SEM;
                        rule_encoding.zscore.pref.ALL.MEAN(:,rule_count) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.ALL.MEAN(:,rule_count) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        
                        rule_encoding.RAW.pref.INC.MEAN(:,rule_count_INC) = temp_block1_MEAN;
                        rule_encoding.RAW.pref.INC.SEM(:,rule_count_INC) = temp_block1_SEM;
                        rule_encoding.RAW.nonpref.INC.MEAN(:,rule_count_INC) = temp_block2_MEAN;
                        rule_encoding.RAW.nonpref.INC.SEM(:,rule_count_INC) = temp_block2_SEM;
                        rule_encoding.zscore.pref.INC.MEAN(:,rule_count_INC) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.INC.MEAN(:,rule_count_INC) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                    elseif mean(FRATE.Cue.Trial_firing_rate_block2) > mean(FRATE.Cue.Trial_firing_rate_block1)
                        rule_encoding.RAW.pref.ALL.MEAN(:,rule_count) = temp_block2_MEAN;
                        rule_encoding.RAW.pref.ALL.SEM(:,rule_count) = temp_block2_SEM;
                        rule_encoding.RAW.nonpref.ALL.MEAN(:,rule_count) = temp_block1_MEAN;
                        rule_encoding.RAW.nonpref.ALL.SEM(:,rule_count) = temp_block1_SEM;
                        rule_encoding.zscore.pref.ALL.MEAN(:,rule_count) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.ALL.MEAN(:,rule_count) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        
                        rule_encoding.RAW.pref.INC.MEAN(:,rule_count_INC) = temp_block2_MEAN;
                        rule_encoding.RAW.pref.INC.SEM(:,rule_count_INC) = temp_block2_SEM;
                        rule_encoding.RAW.nonpref.INC.MEAN(:,rule_count_INC) = temp_block1_MEAN;
                        rule_encoding.RAW.nonpref.INC.SEM(:,rule_count_INC) = temp_block1_SEM;
                        rule_encoding.zscore.pref.INC.MEAN(:,rule_count_INC) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.INC.MEAN(:,rule_count_INC) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                    end
                    rule_count_INC = rule_count_INC + 1;
                elseif mean(FRATE.Task.Trial_firing_rate) < mean(FRATE.Task.Trial_B4_firing_rate)
                    if mean(FRATE.Cue.Trial_firing_rate_block1) > mean(FRATE.Cue.Trial_firing_rate_block2)
                        rule_encoding.RAW.pref.ALL.MEAN(:,rule_count) = temp_block1_MEAN;
                        rule_encoding.RAW.pref.ALL.SEM(:,rule_count) = temp_block1_SEM;
                        rule_encoding.RAW.nonpref.ALL.MEAN(:,rule_count) = temp_block2_MEAN;
                        rule_encoding.RAW.nonpref.ALL.SEM(:,rule_count) = temp_block2_SEM;
                        rule_encoding.zscore.pref.ALL.MEAN(:,rule_count) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.ALL.MEAN(:,rule_count) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        
                        rule_encoding.RAW.pref.DEC.MEAN(:,rule_count_DEC) = temp_block1_MEAN;
                        rule_encoding.RAW.pref.DEC.SEM(:,rule_count_DEC) = temp_block1_SEM;
                        rule_encoding.RAW.nonpref.DEC.MEAN(:,rule_count_DEC) = temp_block2_MEAN;
                        rule_encoding.RAW.nonpref.DEC.SEM(:,rule_count_DEC) = temp_block2_SEM;
                        rule_encoding.zscore.pref.DEC.MEAN(:,rule_count_DEC) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.DEC.MEAN(:,rule_count_DEC) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                    elseif mean(FRATE.Cue.Trial_firing_rate_block2) > mean(FRATE.Cue.Trial_firing_rate_block1)
                        rule_encoding.RAW.pref.ALL.MEAN(:,rule_count) = temp_block2_MEAN;
                        rule_encoding.RAW.pref.ALL.SEM(:,rule_count) = temp_block2_SEM;
                        rule_encoding.RAW.nonpref.ALL.MEAN(:,rule_count) = temp_block1_MEAN;
                        rule_encoding.RAW.nonpref.ALL.SEM(:,rule_count) = temp_block1_SEM;
                        rule_encoding.zscore.pref.ALL.MEAN(:,rule_count) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.ALL.MEAN(:,rule_count) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        
                        rule_encoding.RAW.pref.DEC.MEAN(:,rule_count_DEC) = temp_block2_MEAN;
                        rule_encoding.RAW.pref.DEC.SEM(:,rule_count_DEC) = temp_block2_SEM;
                        rule_encoding.RAW.nonpref.DEC.MEAN(:,rule_count_DEC) = temp_block1_MEAN;
                        rule_encoding.RAW.nonpref.DEC.SEM(:,rule_count_DEC) = temp_block1_SEM;
                        rule_encoding.zscore.pref.DEC.MEAN(:,rule_count_DEC) = (temp_block2_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                        rule_encoding.zscore.nonpref.DEC.MEAN(:,rule_count_DEC) = (temp_block1_MEAN - rule_encoding.zscore.MEAN(rule_count)) / rule_encoding.zscore.STD(rule_count);
                    end
                    rule_count_DEC = rule_count_DEC + 1;
                end
                
                rule_count = rule_count + 1;
                
end

for jj = 1:15001 %generate averaged responses
    if jj <= length(rule_encoding.RAW.pref.ALL.MEAN)
        popSHUFF.rule_encoding.RAW.pref.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.RAW.pref.ALL.MEAN(jj,:));
        popSHUFF.rule_encoding.RAW.pref.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.RAW.pref.ALL.MEAN(jj,:)/sqrt(numel(rule_encoding.RAW.pref.ALL.MEAN(jj,:))-sum(isnan(rule_encoding.RAW.pref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.RAW.nonpref.ALL.MEAN)
        popSHUFF.rule_encoding.RAW.nonpref.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.RAW.nonpref.ALL.MEAN(jj,:));
        popSHUFF.rule_encoding.RAW.nonpref.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.RAW.nonpref.ALL.MEAN(jj,:)/sqrt(numel(rule_encoding.RAW.nonpref.ALL.MEAN(jj,:))-sum(isnan(rule_encoding.RAW.nonpref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.pref.ALL.MEAN)
        popSHUFF.rule_encoding.zscore.pref.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.pref.ALL.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.pref.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.pref.ALL.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.pref.ALL.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.pref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.nonpref.ALL.MEAN)
        popSHUFF.rule_encoding.zscore.nonpref.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.nonpref.ALL.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.nonpref.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.nonpref.ALL.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.nonpref.ALL.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.nonpref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.pref.INC.MEAN)
        popSHUFF.rule_encoding.zscore.pref.INC.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.pref.INC.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.pref.INC.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.pref.INC.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.pref.INC.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.pref.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.nonpref.INC.MEAN)
        popSHUFF.rule_encoding.zscore.nonpref.INC.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.nonpref.INC.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.nonpref.INC.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.nonpref.INC.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.nonpref.INC.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.nonpref.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.pref.DEC.MEAN)
        popSHUFF.rule_encoding.zscore.pref.DEC.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.pref.DEC.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.pref.DEC.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.pref.DEC.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.pref.DEC.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.pref.DEC.MEAN(jj,:)))));
    end
    
    if jj <= length(rule_encoding.zscore.nonpref.DEC.MEAN)
        popSHUFF.rule_encoding.zscore.nonpref.DEC.AVG.MEAN(jj,iShuff) = nanmean(rule_encoding.zscore.nonpref.DEC.MEAN(jj,:));
        popSHUFF.rule_encoding.zscore.nonpref.DEC.AVG.SEM(jj,iShuff) = nanstd(rule_encoding.zscore.nonpref.DEC.MEAN(jj,:)/sqrt(numel(rule_encoding.zscore.nonpref.DEC.MEAN(jj,:))-sum(isnan(rule_encoding.zscore.nonpref.DEC.MEAN(jj,:)))));
    end
end

%% average for location
num_cells = length(find(ALL_matrix(:,2) == 1));
rng('shuffle')
load_vars = datasample(cells_to_use,num_cells,'Replace',false);
% load_vars = randperm(length(ALL_matrix),num_cells);
% load_vars = find(ALL_matrix(:,2) == 1);

location_count = 1;
location_count_INC = 1;
location_count_DEC = 1;

mat_files = dir('*.mat');
for kk = 1:length(load_vars)
    load(mat_files(load_vars(kk)).name);
    mat_overview.fname{kk} = mat_files(load_vars(kk)).name;
%     disp(cat(2,num2str(kk),'/',num2str(length(load_vars))));
        
    location_encoding.zscore.MEAN(location_count) = mean(PETH.Trial.MEAN.trials_PETH);
    location_encoding.zscore.STD(location_count) = std(PETH.Trial.MEAN.trials_PETH);    
    
    arm2_start = find(FRATE.Arm.Trial_firing_rate_groups == 2,1,'first');
    arm3_start = find(FRATE.Arm.Trial_firing_rate_groups == 3,1,'first');
    arm4_start = find(FRATE.Arm.Trial_firing_rate_groups == 4,1,'first');
    
    arm_average_trial_firing_rate(1) = mean(FRATE.Arm.Trial_firing_rate(1:arm2_start-1));
    arm_average_trial_firing_rate(2) = mean(FRATE.Arm.Trial_firing_rate(arm2_start:arm3_start-1));
    arm_average_trial_firing_rate(3) = mean(FRATE.Arm.Trial_firing_rate(arm3_start:arm4_start-1));
    arm_average_trial_firing_rate(4) = mean(FRATE.Arm.Trial_firing_rate(arm4_start:end));
    [B, idx] = sort(arm_average_trial_firing_rate);
    
    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
        case 1
            
            for i = 1:4
                
                switch idx(i)
                    case 1
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor1;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor1;
                    case 2
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor2;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor2;
                    case 3
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor3;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor3;
                    case 4
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor4;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor4;
                end
                
                switch i
                    case 1
                        location_encoding.RAW.one.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.one.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.one.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.one.INC.MEAN(:,location_count_INC) = temp_FR_MEAN;
                        location_encoding.RAW.one.INC.SEM(:,location_count_INC) = temp_FR_SEM;
                        location_encoding.zscore.one.INC.MEAN(:,location_count_INC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 2
                        location_encoding.RAW.two.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.two.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.two.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.two.INC.MEAN(:,location_count_INC) = temp_FR_MEAN;
                        location_encoding.RAW.two.INC.SEM(:,location_count_INC) = temp_FR_SEM;
                        location_encoding.zscore.two.INC.MEAN(:,location_count_INC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 3
                        location_encoding.RAW.three.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.three.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.three.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.three.INC.MEAN(:,location_count_INC) = temp_FR_MEAN;
                        location_encoding.RAW.three.INC.SEM(:,location_count_INC) = temp_FR_SEM;
                        location_encoding.zscore.three.INC.MEAN(:,location_count_INC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 4
                        location_encoding.RAW.four.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.four.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.four.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.four.INC.MEAN(:,location_count_INC) = temp_FR_MEAN;
                        location_encoding.RAW.four.INC.SEM(:,location_count_INC) = temp_FR_SEM;
                        location_encoding.zscore.four.INC.MEAN(:,location_count_INC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                end
                
            end
            location_count_INC = location_count_INC + 1;
            
        case 0
            for i = 1:4
                
                switch idx(i)
                    case 1
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor1;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor1;
                    case 2
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor2;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor2;
                    case 3
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor3;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor3;
                    case 4
                        temp_FR_MEAN = PETH.Arm.MEAN.photosensor4;
                        temp_FR_SEM = PETH.Arm.SEM.photosensor4;
                end
                
                switch i
                    case 1
                        location_encoding.RAW.one.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.one.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.one.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.one.DEC.MEAN(:,location_count_DEC) = temp_FR_MEAN;
                        location_encoding.RAW.one.DEC.SEM(:,location_count_DEC) = temp_FR_SEM;
                        location_encoding.zscore.one.DEC.MEAN(:,location_count_DEC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 2
                        location_encoding.RAW.two.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.two.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.two.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.two.DEC.MEAN(:,location_count_DEC) = temp_FR_MEAN;
                        location_encoding.RAW.two.DEC.SEM(:,location_count_DEC) = temp_FR_SEM;
                        location_encoding.zscore.two.DEC.MEAN(:,location_count_DEC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 3
                        location_encoding.RAW.three.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.three.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.three.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.three.DEC.MEAN(:,location_count_DEC) = temp_FR_MEAN;
                        location_encoding.RAW.three.DEC.SEM(:,location_count_DEC) = temp_FR_SEM;
                        location_encoding.zscore.three.DEC.MEAN(:,location_count_DEC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                    case 4
                        location_encoding.RAW.four.ALL.MEAN(:,location_count) = temp_FR_MEAN;
                        location_encoding.RAW.four.ALL.SEM(:,location_count) = temp_FR_SEM;
                        location_encoding.zscore.four.ALL.MEAN(:,location_count) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                        
                        location_encoding.RAW.four.DEC.MEAN(:,location_count_DEC) = temp_FR_MEAN;
                        location_encoding.RAW.four.DEC.SEM(:,location_count_DEC) = temp_FR_SEM;
                        location_encoding.zscore.four.DEC.MEAN(:,location_count_DEC) = (temp_FR_MEAN - location_encoding.zscore.MEAN(location_count)) / location_encoding.zscore.STD(location_count);
                end
                
            end
            location_count_DEC = location_count_DEC + 1;
            
    end
    location_count = location_count + 1;        
end

for jj = 1:15001 %generate averaged responses
    if jj <= length(location_encoding.RAW.one.ALL.MEAN)
        popSHUFF.location_encoding.RAW.one.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.RAW.one.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.RAW.one.AVG.SEM(jj,iShuff) = nanstd(location_encoding.RAW.one.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.RAW.one.ALL.MEAN(jj,:))-sum(isnan(location_encoding.RAW.one.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.one.ALL.MEAN)
        popSHUFF.location_encoding.zscore.one.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.one.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.one.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.one.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.one.ALL.MEAN(jj,:))-sum(isnan(location_encoding.zscore.one.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.one.INC.MEAN)
        popSHUFF.location_encoding.zscore.one.INC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.one.INC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.one.INC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.one.INC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.one.INC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.one.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.one.DEC.MEAN)
        popSHUFF.location_encoding.zscore.one.DEC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.one.DEC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.one.DEC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.one.DEC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.one.DEC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.one.DEC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.RAW.two.ALL.MEAN)
        popSHUFF.location_encoding.RAW.two.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.RAW.two.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.RAW.two.AVG.SEM(jj,iShuff) = nanstd(location_encoding.RAW.two.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.RAW.two.ALL.MEAN(jj,:))-sum(isnan(location_encoding.RAW.two.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.two.ALL.MEAN)
        popSHUFF.location_encoding.zscore.two.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.two.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.two.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.two.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.two.ALL.MEAN(jj,:))-sum(isnan(location_encoding.zscore.two.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.two.INC.MEAN)
        popSHUFF.location_encoding.zscore.two.INC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.two.INC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.two.INC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.two.INC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.two.INC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.two.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.two.DEC.MEAN)
        popSHUFF.location_encoding.zscore.two.DEC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.two.DEC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.two.DEC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.two.DEC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.two.DEC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.two.DEC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.RAW.three.ALL.MEAN)
        popSHUFF.location_encoding.RAW.three.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.RAW.three.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.RAW.three.AVG.SEM(jj,iShuff) = nanstd(location_encoding.RAW.three.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.RAW.three.ALL.MEAN(jj,:))-sum(isnan(location_encoding.RAW.three.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.three.ALL.MEAN)
        popSHUFF.location_encoding.zscore.three.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.three.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.three.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.three.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.three.ALL.MEAN(jj,:))-sum(isnan(location_encoding.zscore.three.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.three.INC.MEAN)
        popSHUFF.location_encoding.zscore.three.INC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.three.INC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.three.INC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.three.INC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.three.INC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.three.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.three.DEC.MEAN)
        popSHUFF.location_encoding.zscore.three.DEC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.three.DEC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.three.DEC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.three.DEC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.three.DEC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.three.DEC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.RAW.four.ALL.MEAN)
        popSHUFF.location_encoding.RAW.four.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.RAW.four.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.RAW.four.AVG.SEM(jj,iShuff) = nanstd(location_encoding.RAW.four.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.RAW.four.ALL.MEAN(jj,:))-sum(isnan(location_encoding.RAW.four.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.four.ALL.MEAN)
        popSHUFF.location_encoding.zscore.four.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.four.ALL.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.four.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.four.ALL.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.four.ALL.MEAN(jj,:))-sum(isnan(location_encoding.zscore.four.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.four.INC.MEAN)
        popSHUFF.location_encoding.zscore.four.INC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.four.INC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.four.INC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.four.INC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.four.INC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.four.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(location_encoding.zscore.four.DEC.MEAN)
        popSHUFF.location_encoding.zscore.four.DEC.AVG.MEAN(jj,iShuff) = nanmean(location_encoding.zscore.four.DEC.MEAN(jj,:));
        popSHUFF.location_encoding.zscore.four.DEC.AVG.SEM(jj,iShuff) = nanstd(location_encoding.zscore.four.DEC.MEAN(jj,:)/sqrt(numel(location_encoding.zscore.four.DEC.MEAN(jj,:))-sum(isnan(location_encoding.zscore.four.DEC.MEAN(jj,:)))));
    end
    
end

%% outcome
num_cells = length(find(ALL_matrix(:,3) == 1));
rng('shuffle')
load_vars = datasample(cells_to_use,num_cells,'Replace',false);
% load_vars = randperm(length(ALL_matrix),num_cells);
% load_vars = find(ALL_matrix(:,3) == 1);

outcome_count = 1;
outcome_count_INC = 1;
outcome_count_DEC = 1;

mat_files = dir('*.mat');
for kk = 1:length(load_vars)
    load(mat_files(load_vars(kk)).name);
    mat_overview.fname{kk} = mat_files(load_vars(kk)).name;
%     disp(cat(2,num2str(kk),'/',num2str(length(load_vars))));   
                
                outcome_encoding.zscore.MEAN(outcome_count) = mean(PETH.Trial.MEAN.trials_PETH);
                outcome_encoding.zscore.STD(outcome_count) = std(PETH.Trial.MEAN.trials_PETH);
                
                if mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                    if mean(FRATE.Reward.Trial_firing_rate_reward) > mean(FRATE.Reward.Trial_firing_rate_unreward)
                        outcome_encoding.RAW.pref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.pref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.zscore.pref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        
                        outcome_encoding.RAW.pref.INC.MEAN(:,outcome_count_INC) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.pref.INC.SEM(:,outcome_count_INC) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.INC.MEAN(:,outcome_count_INC) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.INC.SEM(:,outcome_count_INC) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.zscore.pref.INC.MEAN(:,outcome_count_INC) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.INC.MEAN(:,outcome_count_INC) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                    elseif mean(FRATE.Reward.Trial_firing_rate_unreward) > mean(FRATE.Reward.Trial_firing_rate_reward)
                        outcome_encoding.RAW.pref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.pref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.zscore.pref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        
                        outcome_encoding.RAW.pref.INC.MEAN(:,outcome_count_INC) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.pref.INC.SEM(:,outcome_count_INC) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.INC.MEAN(:,outcome_count_INC) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.INC.SEM(:,outcome_count_INC) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.zscore.pref.INC.MEAN(:,outcome_count_INC) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.INC.MEAN(:,outcome_count_INC) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                    end
                    outcome_count_INC = outcome_count_INC + 1;
                elseif mean(FRATE.Task.Trial_firing_rate) < mean(FRATE.Task.Trial_B4_firing_rate)
                    if mean(FRATE.Reward.Trial_firing_rate_reward) > mean(FRATE.Reward.Trial_firing_rate_unreward)
                        outcome_encoding.RAW.pref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.pref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.zscore.pref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        
                        outcome_encoding.RAW.pref.DEC.MEAN(:,outcome_count_DEC) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.pref.DEC.SEM(:,outcome_count_DEC) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.DEC.MEAN(:,outcome_count_DEC) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.DEC.SEM(:,outcome_count_DEC) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.zscore.pref.DEC.MEAN(:,outcome_count_DEC) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.DEC.MEAN(:,outcome_count_DEC) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                    elseif mean(FRATE.Reward.Trial_firing_rate_unreward) > mean(FRATE.Reward.Trial_firing_rate_reward)
                        outcome_encoding.RAW.pref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.pref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.MEAN(:,outcome_count) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.ALL.SEM(:,outcome_count) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.zscore.pref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.ALL.MEAN(:,outcome_count) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        
                        outcome_encoding.RAW.pref.DEC.MEAN(:,outcome_count_DEC) = PETH.Trial.MEAN.trials_unrew_PETH;
                        outcome_encoding.RAW.pref.DEC.SEM(:,outcome_count_DEC) = PETH.Trial.SEM.trials_unrew_PETH;
                        outcome_encoding.RAW.nonpref.DEC.MEAN(:,outcome_count_DEC) = PETH.Trial.MEAN.trials_rew_PETH;
                        outcome_encoding.RAW.nonpref.DEC.SEM(:,outcome_count_DEC) = PETH.Trial.SEM.trials_rew_PETH;
                        outcome_encoding.zscore.pref.DEC.MEAN(:,outcome_count_DEC) = (PETH.Trial.MEAN.trials_unrew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                        outcome_encoding.zscore.nonpref.DEC.MEAN(:,outcome_count_DEC) = (PETH.Trial.MEAN.trials_rew_PETH - outcome_encoding.zscore.MEAN(outcome_count)) / outcome_encoding.zscore.STD(outcome_count);
                    end
                    outcome_count_DEC = outcome_count_DEC + 1;
                end
                
                outcome_count = outcome_count + 1;
                
end

for jj = 1:15001 %generate averaged responses
    if jj <= length(outcome_encoding.RAW.pref.ALL.MEAN)
        popSHUFF.outcome_encoding.RAW.pref.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.RAW.pref.ALL.MEAN(jj,:));
        popSHUFF.outcome_encoding.RAW.pref.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.RAW.pref.ALL.MEAN(jj,:)/sqrt(numel(outcome_encoding.RAW.pref.ALL.MEAN(jj,:))-sum(isnan(outcome_encoding.RAW.pref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.RAW.nonpref.ALL.MEAN)
        popSHUFF.outcome_encoding.RAW.nonpref.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.RAW.nonpref.ALL.MEAN(jj,:));
        popSHUFF.outcome_encoding.RAW.nonpref.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.RAW.nonpref.ALL.MEAN(jj,:)/sqrt(numel(outcome_encoding.RAW.nonpref.ALL.MEAN(jj,:))-sum(isnan(outcome_encoding.RAW.nonpref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.pref.ALL.MEAN)
        popSHUFF.outcome_encoding.zscore.pref.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.pref.ALL.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.pref.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.pref.ALL.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.pref.ALL.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.pref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.nonpref.ALL.MEAN)
        popSHUFF.outcome_encoding.zscore.nonpref.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.nonpref.ALL.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.nonpref.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.nonpref.ALL.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.nonpref.ALL.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.nonpref.ALL.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.pref.INC.MEAN)
        popSHUFF.outcome_encoding.zscore.pref.INC.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.pref.INC.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.pref.INC.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.pref.INC.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.pref.INC.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.pref.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.nonpref.INC.MEAN)
        popSHUFF.outcome_encoding.zscore.nonpref.INC.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.nonpref.INC.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.nonpref.INC.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.nonpref.INC.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.nonpref.INC.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.nonpref.INC.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.pref.DEC.MEAN)
        popSHUFF.outcome_encoding.zscore.pref.DEC.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.pref.DEC.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.pref.DEC.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.pref.DEC.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.pref.DEC.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.pref.DEC.MEAN(jj,:)))));
    end
    
    if jj <= length(outcome_encoding.zscore.nonpref.DEC.MEAN)
        popSHUFF.outcome_encoding.zscore.nonpref.DEC.AVG.MEAN(jj,iShuff) = nanmean(outcome_encoding.zscore.nonpref.DEC.MEAN(jj,:));
        popSHUFF.outcome_encoding.zscore.nonpref.DEC.AVG.SEM(jj,iShuff) = nanstd(outcome_encoding.zscore.nonpref.DEC.MEAN(jj,:)/sqrt(numel(outcome_encoding.zscore.nonpref.DEC.MEAN(jj,:))-sum(isnan(outcome_encoding.zscore.nonpref.DEC.MEAN(jj,:)))));
    end
end

%%
clearvars -except cells_to_use popSHUFF ALL_matrix
end

%% plot for modality (normalized)
peak_value = [];
min_value = [];
peak_value(1) = max(rule_encoding.zscore.pref.INC.AVG.MEAN(4001:7000));
peak_value(2) = max(rule_encoding.zscore.nonpref.INC.AVG.MEAN(4001:7000));
maximum_value = max(peak_value);
min_value(1) = min(rule_encoding.zscore.pref.DEC.AVG.MEAN(4001:7000));
min_value(2) = min(rule_encoding.zscore.nonpref.DEC.AVG.MEAN(4001:7000));
minimum_value = min(min_value);

pref_time = -5:.001:10;
pref_time = pref_time(1:length(rule_encoding.zscore.pref.AVG.MEAN));
nonpref_time = -5:.001:10;
nonpref_time = nonpref_time(1:length(rule_encoding.zscore.nonpref.AVG.MEAN));
figure; 
% subplot(3,3,1)
% shadedErrorBar(pref_time,rule_encoding.zscore.pref.AVG.MEAN,rule_encoding.zscore.pref.AVG.SEM,'-r',1);
% hold on; plot(0,-5:.5:50,'.','color','black');
% shadedErrorBar(nonpref_time,rule_encoding.zscore.nonpref.AVG.MEAN,rule_encoding.zscore.nonpref.AVG.SEM,'-b',1);
%   xlim([-5 10]);
% ylim([minimum_value-.5 maximum_value+.5]);
% %     set(gca,'XTick',[]);
% box off;
% xlabel('Time from cue onset (s)');
% ylabel('Normalized firing rate');
% title('All');

subplot(3,2,1)
shadedErrorBar(pref_time,rule_encoding.zscore.pref.INC.AVG.MEAN,rule_encoding.zscore.pref.INC.AVG.SEM,'-r',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(nonpref_time,rule_encoding.zscore.nonpref.INC.AVG.MEAN,rule_encoding.zscore.nonpref.INC.AVG.SEM,'-b',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
xlabel('Time from cue onset (s)');
ylabel('Normalized firing rate');
title('Cue-modulated units that increased post cue-onset');
 set(gca,'FontSize',16);

subplot(3,2,2)
shadedErrorBar(pref_time,rule_encoding.zscore.pref.DEC.AVG.MEAN,rule_encoding.zscore.pref.DEC.AVG.SEM,'-r',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(nonpref_time,rule_encoding.zscore.nonpref.DEC.AVG.MEAN,rule_encoding.zscore.nonpref.DEC.AVG.SEM,'-b',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
 xlabel('Time from cue onset (s)');
% ylabel('Firing rate (Hz)');
ylabel('Normalized firing rate');
%  set(gca,'YTick',[]);
title('Cue-modulated units that decreased post cue-onset');
 set(gca,'FontSize',16);

%% plot for location
peak_value = [];
min_value = [];
peak_value(1) = max(location_encoding.zscore.one.INC.AVG.MEAN(4001:7000));
peak_value(2) = max(location_encoding.zscore.two.INC.AVG.MEAN(4001:7000));
peak_value(3) = max(location_encoding.zscore.three.INC.AVG.MEAN(4001:7000));
peak_value(4) = max(location_encoding.zscore.four.INC.AVG.MEAN(4001:7000));
maximum_value = max(peak_value);
min_value(1) = min(location_encoding.zscore.one.DEC.AVG.MEAN(4001:7000));
min_value(2) = min(location_encoding.zscore.two.DEC.AVG.MEAN(4001:7000));
min_value(3) = min(location_encoding.zscore.three.DEC.AVG.MEAN(4001:7000));
min_value(4) = min(location_encoding.zscore.four.DEC.AVG.MEAN(4001:7000));
minimum_value = min(min_value);

one_time = -5:.001:10;
one_time = one_time(1:length(location_encoding.zscore.one.AVG.MEAN));
two_time = -5:.001:10;
two_time = two_time(1:length(location_encoding.zscore.two.AVG.MEAN));
three_time = -5:.001:10;
three_time = three_time(1:length(location_encoding.zscore.three.AVG.MEAN));
four_time = -5:.001:10;
four_time = four_time(1:length(location_encoding.zscore.four.AVG.MEAN));
% figure; 
% subplot(3,3,4)
% shadedErrorBar(one_time,location_encoding.zscore.one.AVG.MEAN,location_encoding.zscore.one.AVG.SEM,'-m',1);
% hold on; plot(0,-5:.5:50,'.','color','black');
% shadedErrorBar(two_time,location_encoding.zscore.two.AVG.MEAN,location_encoding.zscore.two.AVG.SEM,'-b',1);
% shadedErrorBar(three_time,location_encoding.zscore.three.AVG.MEAN,location_encoding.zscore.three.AVG.SEM,'-k',1);
% shadedErrorBar(four_time,location_encoding.zscore.four.AVG.MEAN,location_encoding.zscore.four.AVG.SEM,'-y',1);
%   xlim([-5 10]);
% ylim([minimum_value-.5 maximum_value+.5]);
% %     set(gca,'XTick',[]);
% box off;
% % xlabel('Time from cue onset (s)');
% ylabel('Normalized firing rate');
% % title('All');

subplot(3,2,3)
shadedErrorBar(one_time,location_encoding.zscore.one.INC.AVG.MEAN,location_encoding.zscore.one.INC.AVG.SEM,'-m',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(two_time,location_encoding.zscore.two.INC.AVG.MEAN,location_encoding.zscore.two.INC.AVG.SEM,'-b',1);
shadedErrorBar(three_time,location_encoding.zscore.three.INC.AVG.MEAN,location_encoding.zscore.three.INC.AVG.SEM,'-k',1);
shadedErrorBar(four_time,location_encoding.zscore.four.INC.AVG.MEAN,location_encoding.zscore.four.INC.AVG.SEM,'-y',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
 xlabel('Time from cue onset (s)');
 ylabel('Normalized firing rate');
% title('Increased');
 set(gca,'FontSize',16);

subplot(3,2,4)
shadedErrorBar(one_time,location_encoding.zscore.one.DEC.AVG.MEAN,location_encoding.zscore.one.DEC.AVG.SEM,'-m',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(two_time,location_encoding.zscore.two.DEC.AVG.MEAN,location_encoding.zscore.two.DEC.AVG.SEM,'-b',1);
shadedErrorBar(three_time,location_encoding.zscore.three.DEC.AVG.MEAN,location_encoding.zscore.three.DEC.AVG.SEM,'-k',1);
shadedErrorBar(four_time,location_encoding.zscore.four.DEC.AVG.MEAN,location_encoding.zscore.four.DEC.AVG.SEM,'-y',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
ylabel('Normalized firing rate');
%  set(gca,'YTick',[]);
 xlabel('Time from cue onset (s)');
% ylabel('Firing rate (Hz)');
% title('Decreased');
 set(gca,'FontSize',16);

%% plot for outcome (normalized)
peak_value = [];
min_value = [];
peak_value(1) = max(outcome_encoding.zscore.pref.INC.AVG.MEAN(4001:7000));
peak_value(2) = max(outcome_encoding.zscore.nonpref.INC.AVG.MEAN(4001:7000));
maximum_value = max(peak_value);
min_value(1) = min(outcome_encoding.zscore.pref.DEC.AVG.MEAN(4001:7000));
min_value(2) = min(outcome_encoding.zscore.nonpref.DEC.AVG.MEAN(4001:7000));
minimum_value = min(min_value);

pref_time = -5:.001:10;
pref_time = pref_time(1:length(outcome_encoding.zscore.pref.AVG.MEAN));
nonpref_time = -5:.001:10;
nonpref_time = nonpref_time(1:length(outcome_encoding.zscore.nonpref.AVG.MEAN));
% figure; 
% subplot(3,3,7)
% shadedErrorBar(pref_time,outcome_encoding.zscore.pref.AVG.MEAN,outcome_encoding.zscore.pref.AVG.SEM,'-r',1);
% hold on; plot(0,-5:.5:50,'.','color','black');
% shadedErrorBar(nonpref_time,outcome_encoding.zscore.nonpref.AVG.MEAN,outcome_encoding.zscore.nonpref.AVG.SEM,'-g',1);
%   xlim([-5 10]);
% ylim([minimum_value-.5 maximum_value+.5]);
% %     set(gca,'XTick',[]);
% box off;
% xlabel('Time from cue onset (s)');
% ylabel('Normalized firing rate');
% % title('All');

subplot(3,2,5)
shadedErrorBar(pref_time,outcome_encoding.zscore.pref.INC.AVG.MEAN,outcome_encoding.zscore.pref.INC.AVG.SEM,'-r',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(nonpref_time,outcome_encoding.zscore.nonpref.INC.AVG.MEAN,outcome_encoding.zscore.nonpref.INC.AVG.SEM,'-g',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
xlabel('Time from cue onset (s)');
ylabel('Normalized firing rate');
% title('Increased');
set(gca,'FontSize',16);

subplot(3,2,6)
shadedErrorBar(pref_time,outcome_encoding.zscore.pref.DEC.AVG.MEAN,outcome_encoding.zscore.pref.DEC.AVG.SEM,'-r',1);
hold on; plot(0,-5:.5:50,'.','color','black');
shadedErrorBar(nonpref_time,outcome_encoding.zscore.nonpref.DEC.AVG.MEAN,outcome_encoding.zscore.nonpref.DEC.AVG.SEM,'-g',1);
  xlim([-5 10]);
ylim([minimum_value-.5 maximum_value+.5]);
%     set(gca,'XTick',[]);
box off;
xlabel('Time from cue onset (s)');
ylabel('Normalized firing rate');
%  set(gca,'YTick',[]);
 set(gca,'FontSize',16);
 
% title('Decreased');