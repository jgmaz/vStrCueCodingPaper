% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

cfg.start = 4001;
cfg.end = 6000;
count = 1;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
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
    
    switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
    
    mat_files_drift{count} = mat_files(kk).name; 
    
    %%
    t = [0 4000];
    binsize = 0.001;
    tbin_edges = t(1):binsize:t(2);
    tbin_centers = tbin_edges(1:end-1)+binsize/2;
    
    spk_count = histc(spk_t,tbin_edges);
    spk_count = spk_count(1:end-1);
    
    binsize = 0.001; % in seconds
    gauss_window = 1./binsize; % 1 second window
    gauss_SD = 0.1./binsize; % 0.1 seconds (100ms) SD
    gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
    gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
    
    temp = -5000;
    photosensor1 = [];
    photosensor1_count = 1;
    photosensor2 = [];
    photosensor2_count = 1;
    
    new_v_old = strcmp(mat_overview.fname{kk}(1:4),'R060');
    switch new_v_old
        case 0
            %% old rats (R053,R056,R057)
            for ik = 1:length(metadata.TrialInfo_block1.trialT)
                if metadata.TrialInfo_block1.photosensorID(ik) == 1%add to count of location of trial
                    for jk = 1:15001
                        if ik == 1
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                                if metadata.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > metadata.TrialInfo_block1.summary(ik-1,9)
                                        photosensor1(jk,photosensor1_count) = NaN;
                                    else
                                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    photosensor1(jk,photosensor1_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo_block1.summary(ik,9)
                                photosensor1(jk:15001,photosensor1_count) = NaN;
                                break
                            else
                                photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor1_count = photosensor1_count + 1;
                elseif metadata.TrialInfo_block1.photosensorID(ik) == 2%add to count of location of trial
                    for jk = 1:15001
                        if ik == 1
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                                if metadata.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > metadata.TrialInfo_block1.summary(ik-1,9)
                                        photosensor2(jk,photosensor2_count) = NaN;
                                    else
                                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    photosensor2(jk,photosensor2_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo_block1.summary(ik,9)
                                photosensor2(jk:15001,photosensor2_count) = NaN;
                                break
                            else
                                photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor2_count = photosensor2_count + 1;
                end
                temp = -5000;
            end
            
            for ip = 1:length(metadata.TrialInfo_block2.trialT)
                if metadata.TrialInfo_block2.photosensorID(ip) == 1%add to count of location of trial
                    for jk = 1:15001
                        if ip == 1
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                                if metadata.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > metadata.TrialInfo_block2.summary(ip-1,9)
                                        photosensor1(jk,photosensor1_count) = NaN;
                                    else
                                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                                    end
                                else
                                    photosensor1(jk,photosensor1_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo_block2.summary(ip,9)
                                photosensor1(jk:15001,photosensor1_count) = NaN;
                                break
                            else
                                photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor1_count = photosensor1_count + 1;
                elseif metadata.TrialInfo_block2.photosensorID(ip) == 2%add to count of location of trial
                    for jk = 1:15001
                        if ip == 1
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                                if metadata.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > metadata.TrialInfo_block2.summary(ip-1,9)
                                        photosensor2(jk,photosensor2_count) = NaN;
                                    else
                                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                                    end
                                else
                                    photosensor2(jk,photosensor2_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo_block2.summary(ip,9)
                                photosensor2(jk:15001,photosensor2_count) = NaN;
                                break
                            else
                                photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo_block1.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor2_count = photosensor2_count + 1;   
                end
                temp = -5000;
            end
            
        case 1 %R060
            for ik = 1:length(metadata.TrialInfo{1}.trialT)
                if metadata.TrialInfo{1}.photosensorID(ik) == 1%add to count of location of trial
                    for jk = 1:15001
                        if ik == 1
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{1}.unnosepoke_to_trialT(ik-1)
                                if metadata.TrialInfo{1}.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > metadata.TrialInfo{1}.summary(ik-1,9)
                                        photosensor1(jk,photosensor1_count) = NaN;
                                    else
                                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    photosensor1(jk,photosensor1_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo{1}.summary(ik,9)
                                photosensor1(jk:15001,photosensor1_count) = NaN;
                                break
                            else
                                photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor1_count = photosensor1_count + 1;
                elseif metadata.TrialInfo{1}.photosensorID(ik) == 2%add to count of location of trial
                    for jk = 1:15001
                        if ik == 1
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{1}.unnosepoke_to_trialT(ik-1)
                                if metadata.TrialInfo{1}.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > metadata.TrialInfo{1}.summary(ik-1,9)
                                        photosensor2(jk,photosensor2_count) = NaN;
                                    else
                                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    photosensor2(jk,photosensor2_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo{1}.summary(ik,9)
                                photosensor2(jk:15001,photosensor2_count) = NaN;
                                break
                            else
                                photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor2_count = photosensor2_count + 1;
                end
                temp = -5000;
            end
            
            for ip = 1:length(metadata.TrialInfo{2}.trialT)
                if metadata.TrialInfo{2}.photosensorID(ip) == 1%add to count of location of trial
                    for jk = 1:15001
                        if ip == 1
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{2}.unnosepoke_to_trialT(ip-1)
                                if metadata.TrialInfo{2}.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > metadata.TrialInfo{2}.summary(ip-1,9)
                                        photosensor1(jk,photosensor1_count) = NaN;
                                    else
                                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                                    end
                                else
                                    photosensor1(jk,photosensor1_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo{2}.summary(ip,9)
                                photosensor1(jk:15001,photosensor1_count) = NaN;
                                break
                            else
                                photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor1_count = photosensor1_count + 1;
                elseif metadata.TrialInfo{2}.photosensorID(ip) == 2%add to count of location of trial
                    for jk = 1:15001
                        if ip == 1
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{2}.unnosepoke_to_trialT(ip-1)
                                if metadata.TrialInfo{2}.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > metadata.TrialInfo{2}.summary(ip-1,9)
                                        photosensor2(jk,photosensor2_count) = NaN;
                                    else
                                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                                    end
                                else
                                    photosensor2(jk,photosensor2_count) = NaN;
                                end
                            elseif temp / 1000 > metadata.TrialInfo{2}.summary(ip,9)
                                photosensor2(jk:15001,photosensor2_count) = NaN;
                                break
                            else
                                photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(metadata.TrialInfo{1}.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    photosensor2_count = photosensor2_count + 1;   
                end
                temp = -5000;
            end
    end
    PETH.Arm.ALL.photosensor1 = photosensor1;
    PETH.Arm.ALL.photosensor2 = photosensor2;
    
    %%    
    rng('shuffle')
    num_trials = round(length(PETH.Trial.ALL.trials_light_PETH(1,:)));
    temp_idx = randperm(num_trials);
    half_1_light = temp_idx(1:num_trials/2);
    half_2_light = temp_idx(num_trials/2+1:end);
    
    num_trials = round(length(PETH.Trial.ALL.trials_sound_PETH(1,:)));
    temp_idx = randperm(num_trials);
    half_1_sound = temp_idx(1:num_trials/2);
    half_2_sound = temp_idx(num_trials/2+1:end);
    
    num_trials = round(length(PETH.Trial.ALL.trials_rew_PETH(1,:)));
    temp_idx = randperm(num_trials);
    half_1_rew = temp_idx(1:num_trials/2);
    half_2_rew = temp_idx(num_trials/2+1:end);
    
    num_trials = round(length(PETH.Trial.ALL.trials_unrew_PETH(1,:)));
    temp_idx = randperm(num_trials);
    half_1_unrew = temp_idx(1:num_trials/2);
    half_2_unrew = temp_idx(num_trials/2+1:end);
    
    num_trials = round(length(PETH.Arm.ALL.photosensor1(1,:)));
    temp_idx = randperm(num_trials);
    half_1_arm1 = temp_idx(1:num_trials/2);
    half_2_arm1 = temp_idx(num_trials/2+1:end);
    
    num_trials = round(length(PETH.Arm.ALL.photosensor2(1,:)));
    temp_idx = randperm(num_trials);
    half_1_arm2 = temp_idx(1:num_trials/2);
    half_2_arm2 = temp_idx(num_trials/2+1:end);
    
    for iBin = 1:15001 %generate averaged responses
%         if iBin <= length(PETH.Trial.ALL.trials_rew_PETH)
    PETHS.mod.trials_light_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_light_PETH(iBin,half_1_light));
    PETHS.mod.trials_light_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_light_PETH(iBin,half_2_light));
    
    PETHS.mod.trials_sound_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_sound_PETH(iBin,half_1_sound));
    PETHS.mod.trials_sound_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_sound_PETH(iBin,half_2_sound));

    PETHS.out.trials_rew_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_rew_PETH(iBin,half_1_rew));
    PETHS.out.trials_rew_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_rew_PETH(iBin,half_2_rew));
    
    PETHS.out.trials_unrew_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_unrew_PETH(iBin,half_1_unrew));
    PETHS.out.trials_unrew_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_unrew_PETH(iBin,half_2_unrew));
    
    PETHS.loc.trials_arm1_1st_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor1(iBin,half_1_arm1));
    PETHS.loc.trials_arm1_2nd_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor1(iBin,half_2_arm1));
    
    PETHS.loc.trials_arm2_1st_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor2(iBin,half_1_arm2));
    PETHS.loc.trials_arm2_2nd_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor2(iBin,half_2_arm2));
%         end
    end
    
    temp = corrcoef(PETHS.mod.trials_light_1st_half{count},PETHS.mod.trials_light_2nd_half{count});
    Corr.mod.ctrl.light(count) = temp(2);
    
    temp = corrcoef(PETHS.mod.trials_sound_1st_half{count},PETHS.mod.trials_sound_2nd_half{count});
    Corr.mod.ctrl.sound(count) = temp(2);
    
    temp = corrcoef(PETHS.mod.trials_light_1st_half{count},PETHS.mod.trials_sound_1st_half{count});
    Corr.mod.cond.l1_v_s1(count) = temp(2);
    
    temp = corrcoef(PETHS.mod.trials_light_1st_half{count},PETHS.mod.trials_sound_2nd_half{count});
    Corr.mod.cond.l1_v_s2(count) = temp(2);
    
    temp = corrcoef(PETHS.mod.trials_light_2nd_half{count},PETHS.mod.trials_sound_1st_half{count});
    Corr.mod.cond.l2_v_s1(count) = temp(2);
    
    temp = corrcoef(PETHS.mod.trials_light_2nd_half{count},PETHS.mod.trials_sound_2nd_half{count});
    Corr.mod.cond.l2_v_s2(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_rew_1st_half{count},PETHS.out.trials_rew_2nd_half{count});
    Corr.out.ctrl.rew(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_unrew_1st_half{count},PETHS.out.trials_unrew_2nd_half{count});
    Corr.out.ctrl.unrew(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_rew_1st_half{count},PETHS.out.trials_unrew_1st_half{count});
    Corr.out.cond.r1_v_un1(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_rew_1st_half{count},PETHS.out.trials_unrew_2nd_half{count});
    Corr.out.cond.r1_v_un2(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_rew_2nd_half{count},PETHS.out.trials_unrew_1st_half{count});
    Corr.out.cond.r2_v_un1(count) = temp(2);
    
    temp = corrcoef(PETHS.out.trials_rew_2nd_half{count},PETHS.out.trials_unrew_2nd_half{count});
    Corr.out.cond.r2_v_un2(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm1_1st_half{count},PETHS.loc.trials_arm1_2nd_half{count});
    Corr.loc.ctrl.arm1(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm2_1st_half{count},PETHS.loc.trials_arm2_2nd_half{count});
    Corr.loc.ctrl.arm2(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm1_1st_half{count},PETHS.loc.trials_arm2_1st_half{count});
    Corr.loc.cond.one1_v_two1(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm1_1st_half{count},PETHS.loc.trials_arm2_2nd_half{count});
    Corr.loc.cond.one1_v_two2(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm1_2nd_half{count},PETHS.loc.trials_arm2_1st_half{count});
    Corr.loc.cond.one2_v_two1(count) = temp(2);
    
    temp = corrcoef(PETHS.loc.trials_arm1_2nd_half{count},PETHS.loc.trials_arm2_2nd_half{count});
    Corr.loc.cond.one2_v_two2(count) = temp(2);
    
    count = count + 1;
    
        case 1
    end
end

%%
[~,Corr.mod.ttest.ctrl.light_sound] = ttest(Corr.mod.ctrl.light,Corr.mod.ctrl.sound);
[~,Corr.mod.ttest.ctrl.cond1_2] = ttest(Corr.mod.cond.l1_v_s1,Corr.mod.cond.l1_v_s2);
[~,Corr.mod.ttest.ctrl.cond2_3] = ttest(Corr.mod.cond.l1_v_s2,Corr.mod.cond.l2_v_s1);
[~,Corr.mod.ttest.ctrl.cond1_4] = ttest(Corr.mod.cond.l1_v_s1,Corr.mod.cond.l2_v_s2);
[~,Corr.mod.ttest.cond.light_cond1] = ttest(Corr.mod.ctrl.light,Corr.mod.cond.l1_v_s1);
[~,Corr.mod.ttest.cond.light_cond2] = ttest(Corr.mod.ctrl.light,Corr.mod.cond.l1_v_s2);
[~,Corr.mod.ttest.cond.light_cond3] = ttest(Corr.mod.ctrl.light,Corr.mod.cond.l2_v_s1);
[~,Corr.mod.ttest.cond.light_cond4] = ttest(Corr.mod.ctrl.light,Corr.mod.cond.l2_v_s2);
[~,Corr.mod.ttest.cond.sound_cond1] = ttest(Corr.mod.ctrl.sound,Corr.mod.cond.l1_v_s1);
[~,Corr.mod.ttest.cond.sound_cond2] = ttest(Corr.mod.ctrl.sound,Corr.mod.cond.l1_v_s2);
[~,Corr.mod.ttest.cond.sound_cond3] = ttest(Corr.mod.ctrl.sound,Corr.mod.cond.l2_v_s1);
[~,Corr.mod.ttest.cond.sound_cond4] = ttest(Corr.mod.ctrl.sound,Corr.mod.cond.l2_v_s2);

[~,Corr.out.ttest.ctrl.rew_unrew] = ttest(Corr.out.ctrl.rew,Corr.out.ctrl.unrew);
[~,Corr.out.ttest.ctrl.cond1_2] = ttest(Corr.out.cond.r1_v_un1,Corr.out.cond.r1_v_un2);
[~,Corr.out.ttest.ctrl.cond2_3] = ttest(Corr.out.cond.r1_v_un2,Corr.out.cond.r2_v_un1);
[~,Corr.out.ttest.ctrl.cond1_4] = ttest(Corr.out.cond.r1_v_un1,Corr.out.cond.r2_v_un2);
[~,Corr.out.ttest.cond.rew_cond1] = ttest(Corr.out.ctrl.rew,Corr.out.cond.r1_v_un1);
[~,Corr.out.ttest.cond.rew_cond2] = ttest(Corr.out.ctrl.rew,Corr.out.cond.r1_v_un2);
[~,Corr.out.ttest.cond.rew_cond3] = ttest(Corr.out.ctrl.rew,Corr.out.cond.r2_v_un1);
[~,Corr.out.ttest.cond.rew_cond4] = ttest(Corr.out.ctrl.rew,Corr.out.cond.r2_v_un2);
[~,Corr.out.ttest.cond.unrew_cond1] = ttest(Corr.out.ctrl.unrew,Corr.out.cond.r1_v_un1);
[~,Corr.out.ttest.cond.unrew_cond2] = ttest(Corr.out.ctrl.unrew,Corr.out.cond.r1_v_un2);
[~,Corr.out.ttest.cond.unrew_cond3] = ttest(Corr.out.ctrl.unrew,Corr.out.cond.r2_v_un1);
[~,Corr.out.ttest.cond.unrew_cond4] = ttest(Corr.out.ctrl.unrew,Corr.out.cond.r2_v_un2);

[~,Corr.loc.ttest.ctrl.arm1_arm2] = ttest(Corr.loc.ctrl.arm1,Corr.loc.ctrl.arm2);
[~,Corr.loc.ttest.ctrl.cond1_2] = ttest(Corr.loc.cond.one1_v_two1,Corr.loc.cond.one1_v_two2);
[~,Corr.loc.ttest.ctrl.cond2_3] = ttest(Corr.loc.cond.one1_v_two2,Corr.loc.cond.one2_v_two1);
[~,Corr.loc.ttest.ctrl.cond1_4] = ttest(Corr.loc.cond.one1_v_two1,Corr.loc.cond.one2_v_two2);
[~,Corr.loc.ttest.cond.arm1_cond1] = ttest(Corr.loc.ctrl.arm1,Corr.loc.cond.one1_v_two1);
[~,Corr.loc.ttest.cond.arm1_cond2] = ttest(Corr.loc.ctrl.arm1,Corr.loc.cond.one1_v_two2);
[~,Corr.loc.ttest.cond.arm1_cond3] = ttest(Corr.loc.ctrl.arm1,Corr.loc.cond.one2_v_two1);
[~,Corr.loc.ttest.cond.arm1_cond4] = ttest(Corr.loc.ctrl.arm1,Corr.loc.cond.one2_v_two2);
[~,Corr.loc.ttest.cond.arm2_cond1] = ttest(Corr.loc.ctrl.arm2,Corr.loc.cond.one1_v_two1);
[~,Corr.loc.ttest.cond.arm2_cond2] = ttest(Corr.loc.ctrl.arm2,Corr.loc.cond.one1_v_two2);
[~,Corr.loc.ttest.cond.arm2_cond3] = ttest(Corr.loc.ctrl.arm2,Corr.loc.cond.one2_v_two1);
[~,Corr.loc.ttest.cond.arm2_cond4] = ttest(Corr.loc.ctrl.arm2,Corr.loc.cond.one2_v_two2);

%%
figure
subplot(3,4,1)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.ctrl.sound))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v ctrl sound; p = ',num2str(Corr.mod.ttest.ctrl.light_sound)))
xlabel('ctrl light'); ylabel('ctrl sound');
subplot(3,4,2)
scatter(abs(Corr.mod.cond.l1_v_s1),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp2; p = ',num2str(Corr.mod.ttest.ctrl.cond1_2)))
xlabel('comparison 1'); ylabel('comparison 2');
subplot(3,4,3)
scatter(abs(Corr.mod.cond.l1_v_s2),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'comp2 v comp3; p = ',num2str(Corr.mod.ttest.ctrl.cond2_3)))
xlabel('comparison 2'); ylabel('comparison 3');
subplot(3,4,4)
scatter(abs(Corr.mod.cond.l1_v_s1),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp4; p = ',num2str(Corr.mod.ttest.ctrl.cond1_4)))
xlabel('comparison 1'); ylabel('comparison 4');

subplot(3,4,5)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l1_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp1; p = ',num2str(Corr.mod.ttest.cond.light_cond1)))
xlabel('ctrl light'); ylabel('comparison 1');
subplot(3,4,6)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp2; p = ',num2str(Corr.mod.ttest.cond.light_cond2)))
xlabel('ctrl light'); ylabel('comparison 2');
subplot(3,4,7)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp3; p = ',num2str(Corr.mod.ttest.cond.light_cond3)))
xlabel('ctrl light'); ylabel('comparison 3');
subplot(3,4,8)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp4; p = ',num2str(Corr.mod.ttest.cond.light_cond4)))
xlabel('ctrl light'); ylabel('comparison 4');

subplot(3,4,9)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l1_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp1; p = ',num2str(Corr.mod.ttest.cond.sound_cond1)))
xlabel('ctrl sound'); ylabel('comparison 1');
subplot(3,4,10)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp2; p = ',num2str(Corr.mod.ttest.cond.sound_cond2)))
xlabel('ctrl sound'); ylabel('comparison 2');
subplot(3,4,11)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp3; p = ',num2str(Corr.mod.ttest.cond.sound_cond3)))
xlabel('ctrl sound'); ylabel('comparison 3');
subplot(3,4,12)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp4; p = ',num2str(Corr.mod.ttest.cond.sound_cond4)))
xlabel('ctrl sound'); ylabel('comparison 4');

%%
figure
subplot(3,4,1)
scatter(abs(Corr.out.ctrl.rew),abs(Corr.out.ctrl.unrew))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl rew v ctrl unrew; p = ',num2str(Corr.out.ttest.ctrl.rew_unrew)))
xlabel('ctrl rew'); ylabel('ctrl unrew');
subplot(3,4,2)
scatter(abs(Corr.out.cond.r1_v_un1),abs(Corr.out.cond.r1_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp2; p = ',num2str(Corr.out.ttest.ctrl.cond1_2)))
xlabel('comparison 1'); ylabel('comparison 2');
subplot(3,4,3)
scatter(abs(Corr.out.cond.r1_v_un2),abs(Corr.out.cond.r2_v_un1))
hold on; line([0 1],[0 1]); title(cat(2,'comp2 v comp3; p = ',num2str(Corr.out.ttest.ctrl.cond2_3)))
xlabel('comparison 2'); ylabel('comparison 3');
subplot(3,4,4)
scatter(abs(Corr.out.cond.r1_v_un1),abs(Corr.out.cond.r2_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp4; p = ',num2str(Corr.out.ttest.ctrl.cond1_4)))
xlabel('comparison 1'); ylabel('comparison 4');

subplot(3,4,5)
scatter(abs(Corr.out.ctrl.rew),abs(Corr.out.cond.r1_v_un1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl rew v comp1; p = ',num2str(Corr.out.ttest.cond.rew_cond1)))
xlabel('ctrl rew'); ylabel('comparison 1');
subplot(3,4,6)
scatter(abs(Corr.out.ctrl.rew),abs(Corr.out.cond.r1_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl rew v comp2; p = ',num2str(Corr.out.ttest.cond.rew_cond2)))
xlabel('ctrl rew'); ylabel('comparison 2');
subplot(3,4,7)
scatter(abs(Corr.out.ctrl.rew),abs(Corr.out.cond.r2_v_un1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl rew v comp3; p = ',num2str(Corr.out.ttest.cond.rew_cond3)))
xlabel('ctrl rew'); ylabel('comparison 3');
subplot(3,4,8)
scatter(abs(Corr.out.ctrl.rew),abs(Corr.out.cond.r2_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl rew v comp4; p = ',num2str(Corr.out.ttest.cond.rew_cond4)))
xlabel('ctrl rew'); ylabel('comparison 4');

subplot(3,4,9)
scatter(abs(Corr.out.ctrl.unrew),abs(Corr.out.cond.r1_v_un1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl unrew v comp1; p = ',num2str(Corr.out.ttest.cond.unrew_cond1)))
xlabel('ctrl unrew'); ylabel('comparison 1');
subplot(3,4,10)
scatter(abs(Corr.out.ctrl.unrew),abs(Corr.out.cond.r1_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl unrew v comp2; p = ',num2str(Corr.out.ttest.cond.unrew_cond2)))
xlabel('ctrl unrew'); ylabel('comparison 2');
subplot(3,4,11)
scatter(abs(Corr.out.ctrl.unrew),abs(Corr.out.cond.r2_v_un1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl unrew v comp3; p = ',num2str(Corr.out.ttest.cond.unrew_cond3)))
xlabel('ctrl unrew'); ylabel('comparison 3');
subplot(3,4,12)
scatter(abs(Corr.out.ctrl.unrew),abs(Corr.out.cond.r2_v_un2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl unrew v comp4; p = ',num2str(Corr.out.ttest.cond.unrew_cond4)))
xlabel('ctrl unrew'); ylabel('comparison 4');

%%
figure
subplot(3,4,1)
scatter(abs(Corr.loc.ctrl.arm1),abs(Corr.loc.ctrl.arm2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm1 v ctrl arm2; p = ',num2str(Corr.loc.ttest.ctrl.arm1_arm2)))
xlabel('ctrl arm1'); ylabel('ctrl arm2');
subplot(3,4,2)
scatter(abs(Corr.loc.cond.one1_v_two1),abs(Corr.loc.cond.one1_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp2; p = ',num2str(Corr.loc.ttest.ctrl.cond1_2)))
xlabel('comparison 1'); ylabel('comparison 2');
subplot(3,4,3)
scatter(abs(Corr.loc.cond.one1_v_two2),abs(Corr.loc.cond.one2_v_two1))
hold on; line([0 1],[0 1]); title(cat(2,'comp2 v comp3; p = ',num2str(Corr.loc.ttest.ctrl.cond2_3)))
xlabel('comparison 2'); ylabel('comparison 3');
subplot(3,4,4)
scatter(abs(Corr.loc.cond.one1_v_two1),abs(Corr.loc.cond.one2_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp4; p = ',num2str(Corr.loc.ttest.ctrl.cond1_4)))
xlabel('comparison 1'); ylabel('comparison 4');

subplot(3,4,5)
scatter(abs(Corr.loc.ctrl.arm1),abs(Corr.loc.cond.one1_v_two1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm1 v comp1; p = ',num2str(Corr.loc.ttest.cond.arm1_cond1)))
xlabel('ctrl arm1'); ylabel('comparison 1');
subplot(3,4,6)
scatter(abs(Corr.loc.ctrl.arm1),abs(Corr.loc.cond.one1_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm1 v comp2; p = ',num2str(Corr.loc.ttest.cond.arm1_cond2)))
xlabel('ctrl arm1'); ylabel('comparison 2');
subplot(3,4,7)
scatter(abs(Corr.loc.ctrl.arm1),abs(Corr.loc.cond.one2_v_two1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm1 v comp3; p = ',num2str(Corr.loc.ttest.cond.arm1_cond3)))
xlabel('ctrl arm1'); ylabel('comparison 3');
subplot(3,4,8)
scatter(abs(Corr.loc.ctrl.arm1),abs(Corr.loc.cond.one2_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm1 v comp4; p = ',num2str(Corr.loc.ttest.cond.arm1_cond4)))
xlabel('ctrl arm1'); ylabel('comparison 4');

subplot(3,4,9)
scatter(abs(Corr.loc.ctrl.arm2),abs(Corr.loc.cond.one1_v_two1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm2 v comp1; p = ',num2str(Corr.loc.ttest.cond.arm2_cond1)))
xlabel('ctrl arm2'); ylabel('comparison 1');
subplot(3,4,10)
scatter(abs(Corr.loc.ctrl.arm2),abs(Corr.loc.cond.one1_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm2 v comp2; p = ',num2str(Corr.loc.ttest.cond.arm2_cond2)))
xlabel('ctrl arm2'); ylabel('comparison 2');
subplot(3,4,11)
scatter(abs(Corr.loc.ctrl.arm2),abs(Corr.loc.cond.one2_v_two1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm2 v comp3; p = ',num2str(Corr.loc.ttest.cond.arm2_cond3)))
xlabel('ctrl arm2'); ylabel('comparison 3');
subplot(3,4,12)
scatter(abs(Corr.loc.ctrl.arm2),abs(Corr.loc.cond.one2_v_two2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl arm2 v comp4; p = ',num2str(Corr.loc.ttest.cond.arm2_cond4)))
xlabel('ctrl arm2'); ylabel('comparison 4');

%% set up variables for LME (mod)
temp_Corr = 1:length(Corr.mod.ctrl.light);
Corr_mod_summary(:,1) = cat(2,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr)';

temp_length = length(Corr.mod.ctrl.light);
for iCell = 1:temp_length
temp_Corr_cells{iCell,1} = 'A';
end
for iCell = temp_length+1:temp_length*2
temp_Corr_cells{iCell,1} = 'B';
end
for iCell = temp_length*2+1:temp_length*3
temp_Corr_cells{iCell,1} = 'C';
end
for iCell = temp_length*3+1:temp_length*4
temp_Corr_cells{iCell,1} = 'D';
end
for iCell = temp_length*4+1:temp_length*5
temp_Corr_cells{iCell,1} = 'E';
end
for iCell = temp_length*5+1:temp_length*6
temp_Corr_cells{iCell,1} = 'F';
end

% temp_Corr_cells(1:temp_length,1) = 'A';
% temp_Corr_cells(temp_length+1:temp_length*2,1) = 'B';
% temp_Corr_cells(temp_length*2+1:temp_length*3,1) = 'C';
% temp_Corr_cells(temp_length*3+1:temp_length*4,1) = 'D';
% temp_Corr_cells(temp_length*4+1:temp_length*5,1) = 'E';
% temp_Corr_cells(temp_length*5+1:temp_length*6,1) = 'F';
% Corr_mod_summary(:,2) = temp_Corr_cells;

% temp_Corr_cond(1:temp_length,1) = 1;
% temp_Corr_cond(temp_length+1:temp_length*2,1) = 1;
% temp_Corr_cond(temp_length*2+1:temp_length*3,1) = 2;
% temp_Corr_cond(temp_length*3+1:temp_length*4,1) = 2;
% temp_Corr_cond(temp_length*4+1:temp_length*5,1) = 2;
% temp_Corr_cond(temp_length*5+1:temp_length*6,1) = 2;
% Corr_mod_summary(:,3) = temp_Corr_cond;

Corr_mod_summary(:,2) = cat(2,Corr.mod.ctrl.light,Corr.mod.ctrl.sound,Corr.mod.cond.l1_v_s1,Corr.mod.cond.l1_v_s2,Corr.mod.cond.l2_v_s1,Corr.mod.cond.l2_v_s2)';

%% linear mixed effects model (mod)
Corr.mod.LME.tbl = table(Corr_mod_summary(:,1),temp_Corr_cells,Corr_mod_summary(:,2),'VariableNames',{'CellID','compType','corrCoeff'});
Corr.mod.LME.tbl.compType = nominal(Corr.mod.LME.tbl.compType);
Corr.mod.LME.lme = fitlme(Corr.mod.LME.tbl,'corrCoeff~compType+(1|CellID)');%,'DummyVarCoding','effects');% +(CueType-1|RatID)'); <- matlab includes this, not sure why %lme for Trial length with cue type as fixed effect and a random intercept for each rat.
Corr.mod.LME.lme_reduced = fitlme(Corr.mod.LME.tbl,'corrCoeff~1+(1|CellID)');

Corr.mod.LME.comparison = compare(Corr.mod.LME.lme_reduced,Corr.mod.LME.lme);

%% set up variables for LME (loc)
temp_Corr = 1:length(Corr.loc.ctrl.light);
Corr_loc_summary(:,1) = cat(2,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr)';

temp_length = length(Corr.loc.ctrl.light);
for iCell = 1:temp_length
temp_Corr_cells{iCell,1} = 'A';
end
for iCell = temp_length+1:temp_length*2
temp_Corr_cells{iCell,1} = 'B';
end
for iCell = temp_length*2+1:temp_length*3
temp_Corr_cells{iCell,1} = 'C';
end
for iCell = temp_length*3+1:temp_length*4
temp_Corr_cells{iCell,1} = 'D';
end
for iCell = temp_length*4+1:temp_length*5
temp_Corr_cells{iCell,1} = 'E';
end
for iCell = temp_length*5+1:temp_length*6
temp_Corr_cells{iCell,1} = 'F';
end

% temp_Corr_cells(1:temp_length,1) = 'A';
% temp_Corr_cells(temp_length+1:temp_length*2,1) = 'B';
% temp_Corr_cells(temp_length*2+1:temp_length*3,1) = 'C';
% temp_Corr_cells(temp_length*3+1:temp_length*4,1) = 'D';
% temp_Corr_cells(temp_length*4+1:temp_length*5,1) = 'E';
% temp_Corr_cells(temp_length*5+1:temp_length*6,1) = 'F';
% Corr_loc_summary(:,2) = temp_Corr_cells;

% temp_Corr_cond(1:temp_length,1) = 1;
% temp_Corr_cond(temp_length+1:temp_length*2,1) = 1;
% temp_Corr_cond(temp_length*2+1:temp_length*3,1) = 2;
% temp_Corr_cond(temp_length*3+1:temp_length*4,1) = 2;
% temp_Corr_cond(temp_length*4+1:temp_length*5,1) = 2;
% temp_Corr_cond(temp_length*5+1:temp_length*6,1) = 2;
% Corr_loc_summary(:,3) = temp_Corr_cond;

Corr_loc_summary(:,2) = cat(2,Corr.loc.ctrl.arm1,Corr.loc.ctrl.arm2,Corr.loc.cond.one1_v_two1,Corr.loc.cond.one1_v_two2,Corr.loc.cond.one2_v_two1,Corr.loc.cond.one2_v_two2)';

%% linear mixed effects model (loc)
Corr.loc.LME.tbl = table(Corr_loc_summary(:,1),temp_Corr_cells,Corr_loc_summary(:,2),'VariableNames',{'CellID','compType','corrCoeff'});
Corr.loc.LME.tbl.compType = nominal(Corr.loc.LME.tbl.compType);
Corr.loc.LME.lme = fitlme(Corr.loc.LME.tbl,'corrCoeff~compType+(1|CellID)');%,'DummyVarCoding','effects');% +(CueType-1|RatID)'); <- matlab includes this, not sure why %lme for Trial length with cue type as fixed effect and a random intercept for each rat.
Corr.loc.LME.lme_reduced = fitlme(Corr.loc.LME.tbl,'corrCoeff~1+(1|CellID)');

Corr.loc.LME.comparison = compare(Corr.loc.LME.lme_reduced,Corr.loc.LME.lme);

%% set up variables for LME (out)
temp_Corr = 1:length(Corr.out.ctrl.light);
Corr_out_summary(:,1) = cat(2,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr)';

temp_length = length(Corr.out.ctrl.light);
for iCell = 1:temp_length
temp_Corr_cells{iCell,1} = 'A';
end
for iCell = temp_length+1:temp_length*2
temp_Corr_cells{iCell,1} = 'B';
end
for iCell = temp_length*2+1:temp_length*3
temp_Corr_cells{iCell,1} = 'C';
end
for iCell = temp_length*3+1:temp_length*4
temp_Corr_cells{iCell,1} = 'D';
end
for iCell = temp_length*4+1:temp_length*5
temp_Corr_cells{iCell,1} = 'E';
end
for iCell = temp_length*5+1:temp_length*6
temp_Corr_cells{iCell,1} = 'F';
end

% temp_Corr_cells(1:temp_length,1) = 'A';
% temp_Corr_cells(temp_length+1:temp_length*2,1) = 'B';
% temp_Corr_cells(temp_length*2+1:temp_length*3,1) = 'C';
% temp_Corr_cells(temp_length*3+1:temp_length*4,1) = 'D';
% temp_Corr_cells(temp_length*4+1:temp_length*5,1) = 'E';
% temp_Corr_cells(temp_length*5+1:temp_length*6,1) = 'F';
% Corr_out_summary(:,2) = temp_Corr_cells;

% temp_Corr_cond(1:temp_length,1) = 1;
% temp_Corr_cond(temp_length+1:temp_length*2,1) = 1;
% temp_Corr_cond(temp_length*2+1:temp_length*3,1) = 2;
% temp_Corr_cond(temp_length*3+1:temp_length*4,1) = 2;
% temp_Corr_cond(temp_length*4+1:temp_length*5,1) = 2;
% temp_Corr_cond(temp_length*5+1:temp_length*6,1) = 2;
% Corr_out_summary(:,3) = temp_Corr_cond;

Corr_out_summary(:,2) = cat(2,Corr.out.ctrl.rew,Corr.out.ctrl.unrew,Corr.out.cond.r1_v_un1,Corr.out.cond.r1_v_un2,Corr.out.cond.r2_v_un1,Corr.out.cond.r2_v_un2)';

%% linear mixed effects model (out)
Corr.out.LME.tbl = table(Corr_out_summary(:,1),temp_Corr_cells,Corr_out_summary(:,2),'VariableNames',{'CellID','compType','corrCoeff'});
Corr.out.LME.tbl.compType = nominal(Corr.out.LME.tbl.compType);
Corr.out.LME.lme = fitlme(Corr.out.LME.tbl,'corrCoeff~compType+(1|CellID)');%,'DummyVarCoding','effects');% +(CueType-1|RatID)'); <- matlab includes this, not sure why %lme for Trial length with cue type as fixed effect and a random intercept for each rat.
Corr.out.LME.lme_reduced = fitlme(Corr.out.LME.tbl,'corrCoeff~1+(1|CellID)');

Corr.out.LME.comparison = compare(Corr.out.LME.lme_reduced,Corr.out.LME.lme);

%% means for table
Corr.mod.mean.light = nanmean(abs(Corr.mod.ctrl.light));
Corr.mod.mean.sound = nanmean(abs(Corr.mod.ctrl.sound));
Corr.mod.mean.cond1 = nanmean(abs(Corr.mod.cond.l1_v_s1));
Corr.mod.mean.cond2 = nanmean(abs(Corr.mod.cond.l1_v_s2));
Corr.mod.mean.cond3 = nanmean(abs(Corr.mod.cond.l2_v_s1));
Corr.mod.mean.cond4 = nanmean(abs(Corr.mod.cond.l2_v_s2));

Corr.loc.mean.arm1 = nanmean(abs(Corr.loc.ctrl.arm1));
Corr.loc.mean.arm2 = nanmean(abs(Corr.loc.ctrl.arm2));
Corr.loc.mean.cond1 = nanmean(abs(Corr.loc.cond.one1_v_two1));
Corr.loc.mean.cond2 = nanmean(abs(Corr.loc.cond.one1_v_two2));
Corr.loc.mean.cond3 = nanmean(abs(Corr.loc.cond.one2_v_two1));
Corr.loc.mean.cond4 = nanmean(abs(Corr.loc.cond.one2_v_two2));

Corr.out.mean.rew = nanmean(abs(Corr.out.ctrl.rew));
Corr.out.mean.unrew = nanmean(abs(Corr.out.ctrl.unrew));
Corr.out.mean.cond1 = nanmean(abs(Corr.out.cond.r1_v_un1));
Corr.out.mean.cond2 = nanmean(abs(Corr.out.cond.r1_v_un2));
Corr.out.mean.cond3 = nanmean(abs(Corr.out.cond.r2_v_un1));
Corr.out.mean.cond4 = nanmean(abs(Corr.out.cond.r2_v_un2));

%% RM ANOVA
%table coeff, 1v2, lvs, subj)
Predictor = {'mod' 'loc' 'out'};
meas = {'A' 'B' 'D' 'E'}; % no C and F as they are 1v1 and 2v2, respectively.

for iPred = 1:3
    for iMeas = 1:length(meas)
        startT = find(Corr.(Predictor{iPred}).LME.tbl.compType == meas{iMeas},1,'first');
        endT = find(Corr.(Predictor{iPred}).LME.tbl.compType == meas{iMeas},1,'last');
        count = 1;
        for iData = startT:endT
            RMA.(Predictor{iPred}).RAW(count,iMeas) = Corr.(Predictor{iPred}).LME.tbl.corrCoeff(iData);
            count = count + 1;
        end
    end
    RMA.(Predictor{iPred}).tbl = table(RMA.(Predictor{iPred}).RAW(:,1),RMA.(Predictor{iPred}).RAW(:,2),RMA.(Predictor{iPred}).RAW(:,3),RMA.(Predictor{iPred}).RAW(:,4),'VariableNames',{'w1' 'w2' 'a1' 'a2'});
end

%%
f1 = [1 2 1 2];
f2 = [1 2 2 1];

for iPred = 1:3
    RMA.(Predictor{iPred}).RAWcat = [];
    RMA.(Predictor{iPred}).RAWf1 = [];
    RMA.(Predictor{iPred}).RAWf2 = [];
    RMA.(Predictor{iPred}).RAWid = [];
    for iMeas = 1:length(meas)
        temp1 = [];
        temp2 = [];
        temp3 = [];
        RMA.(Predictor{iPred}).RAWcat = cat(1,RMA.(Predictor{iPred}).RAWcat,RMA.(Predictor{iPred}).RAW(:,iMeas));
        temp1(1:length(RMA.(Predictor{iPred}).RAW(:,iMeas)),1) = f1(iMeas);
        temp2(1:length(RMA.(Predictor{iPred}).RAW(:,iMeas)),1) = f2(iMeas);
        temp3(:,1) = 1:length(RMA.(Predictor{iPred}).RAW(:,iMeas));
        RMA.(Predictor{iPred}).RAWf1 = cat(1,RMA.(Predictor{iPred}).RAWf1,temp1);
        RMA.(Predictor{iPred}).RAWf2 = cat(1,RMA.(Predictor{iPred}).RAWf2,temp2);
        RMA.(Predictor{iPred}).RAWid = cat(1,RMA.(Predictor{iPred}).RAWid,temp3);       
    end
    RMA.(Predictor{iPred}).RMtbl = table(RMA.(Predictor{iPred}).RAWcat,categorical(RMA.(Predictor{iPred}).RAWf1),categorical(RMA.(Predictor{iPred}).RAWf2),RMA.(Predictor{iPred}).RAWid,...
        'VariableNames',{'Corr' 'factor1' 'factor2' 'cellid'});
%     RMA.(Predictor{iPred}).RM = fitrm(RMA.(Predictor{iPred}).RMtbl,'Corr ~ factor1*factor2','WithinDesign',within)
end

%%
for iPred = 1:3
    [RMA.(Predictor{iPred}).RM.p RMA.(Predictor{iPred}).RM.tbl RMA.(Predictor{iPred}).RM.stats] = anovan(RMA.(Predictor{iPred}).RAWcat,{categorical(RMA.(Predictor{iPred}).RAWf1),categorical(RMA.(Predictor{iPred}).RAWf2),categorical(RMA.(Predictor{iPred}).RAWid)},...
    'model','interaction','random',[3]);
end