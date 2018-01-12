function FRATE = genFRATE(meta,spk_t,dataPoint)
% function FRATE = genFRATE(meta,spk_t,dataPoint)
%
% Generates average, and trial-by-trial firing rates broken down a variety
% of ways
%
% INPUTS:
% meta: meta file containing AMPX and task data
% spk_t: vector of spike times
% dataPoint: list of trial and nosepoke start times (in ~ms)
%
% OUTPUTS:
% FRATE: struct containing pre- and post-trial firing rates

%% shuffle data
tic
disp('bootstrapping')

spk_t_pre = spk_t(spk_t < meta.data_pre.tvec(end)); % first I need to separate the spikes by the session they occured in so I do not get shuffled spikes occuring during intersession periods / keep it within a session
spk_t_block1 = spk_t(spk_t > (meta.data_pre.tvec(end) + meta.gap1.tvec(end)) & spk_t < (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end)));
spk_t_block2 = spk_t(spk_t > (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end)) & spk_t < (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.data_block2.tvec(end)));
spk_t_post = spk_t(spk_t > (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.data_block2.tvec(end) + meta.gap3.tvec(end)));

ISI_pre = diff(spk_t_pre);
ISI_block1 = diff(spk_t_block1);
ISI_block2 = diff(spk_t_block2);
ISI_post = diff(spk_t_post);

for kk = 1:1000
    idx_pre = randperm(length(ISI_pre));
    idx_block1 = randperm(length(ISI_block1));
    idx_block2 = randperm(length(ISI_block2));
    idx_post = randperm(length(ISI_post));
    
    shuff_ISI_pre = ISI_pre(idx_pre);
    shuff_ISI_block1 = ISI_block1(idx_block1);
    shuff_ISI_block2 = ISI_block2(idx_block2);
    shuff_ISI_post = ISI_post(idx_post);
    
    if isempty(shuff_ISI_pre)
        shuff_spk_t_pre = [];
    else
    shuff_spk_t_pre(:,kk) = cumsum(shuff_ISI_pre); % is 1 spike short of original data, doesn't factor in initial time for first spike
    end
    shuff_spk_t_block1(:,kk) = cumsum(shuff_ISI_block1);
    shuff_spk_t_block2(:,kk) = cumsum(shuff_ISI_block2);
    if isempty(shuff_ISI_post)
        shuff_spk_t_post = [];
    else
    shuff_spk_t_post(:,kk) = cumsum(shuff_ISI_post);
    end
end

shuff_spk_t_block1_adjusted = shuff_spk_t_block1 + (meta.data_pre.tvec(end) + meta.gap1.tvec(end));
shuff_spk_t_block2_adjusted = shuff_spk_t_block2 + (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end));
shuff_spk_t_post_adjusted = shuff_spk_t_post + (meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.data_block2.tvec(end) + meta.gap3.tvec(end));

shuff_spk_t = cat(1,shuff_spk_t_pre,shuff_spk_t_block1_adjusted,shuff_spk_t_block2_adjusted,shuff_spk_t_post_adjusted);

%% average firing rate for each block
block1_start = meta.data_pre.tvec(end) + meta.gap1.tvec(end);
block1_end = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end);
block2_start = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end);
block2_end = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.data_block2.tvec(end);

block1_spk = spk_t(spk_t > block1_start & spk_t < block1_end);
block2_spk = spk_t(spk_t > block2_start & spk_t < block2_end);

FRATE.Overall.firing_rate_block1 = length(block1_spk) / meta.data_block1.tvec(end);
FRATE.Overall.firing_rate_block2 = length(block2_spk) / meta.data_block2.tvec(end);
FRATE.Overall.firing_rate_total = (length(block1_spk) + length(block2_spk)) / (meta.data_block1.tvec(end) + meta.data_block2.tvec(end));

%% trials
end_time = nan(length(meta.TrialInfo_block1.trialT),1);
end_time2 = nan(length(meta.TrialInfo_block2.trialT),1);
firing_rate = cat(1,end_time,end_time2);
firing_rate_block1 = nan(length(meta.TrialInfo_block1.trialT),1);
firing_rate_block2 = nan(length(meta.TrialInfo_block2.trialT),1);

for ik = 1:length(meta.TrialInfo_block1.trialT)
    if meta.TrialInfo_block1.trial_length_analysis(ik) < 1
        trial_block1.trial_length(ik) = meta.TrialInfo_block1.trial_length_analysis(ik);
        start_time_trials(ik) = dataPoint.Trials(ik)/1000 - meta.TrialInfo_block1.trial_length_analysis(ik);
    else
        trial_block1.trial_length(ik) = 1;
        start_time_trials(ik) = dataPoint.Trials(ik)/1000 - 1;
    end
    end_time(ik) = dataPoint.Trials(ik)/1000 + trial_block1.trial_length(ik);
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > dataPoint.Trials(ik)/1000 & spk_t < end_time(ik));
    for kl = 1:1000
        these_spk_shuff{kl} = (shuff_spk_t(:,kl) > dataPoint.Trials(ik)/1000 & shuff_spk_t(:,kl) < end_time(ik));
    end
    
    % convert to firing rate and store
    firing_rate(ik) = length(these_spk) / trial_block1.trial_length(ik);
    firing_rate_block1(ik) = length(these_spk) / trial_block1.trial_length(ik);
    for kl = 1:1000
        shuff_firing_rate(ik,kl) = length(these_spk_shuff{kl}) / trial_block1.trial_length(ik);
    end
    
end

for ik = 1:length(meta.TrialInfo_block2.trialT) %-1 % -1 just to do trial minus last one as it isn't always completed
    if meta.TrialInfo_block2.trial_length_analysis(ik) < 1
        trial_block2.trial_length(ik) = meta.TrialInfo_block2.trial_length_analysis(ik);
        start_time_trials(length(meta.TrialInfo_block1.trialT) + ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 - meta.TrialInfo_block2.trial_length_analysis(ik);
    else
        trial_block2.trial_length(ik) = 1;
        start_time_trials(length(meta.TrialInfo_block1.trialT) + ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 - 1;
    end
    end_time2(ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 + trial_block2.trial_length(ik);
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 & spk_t < end_time2(ik));
    for kl = 1:1000
    these_spk_shuff{kl} = shuff_spk_t(shuff_spk_t(:,kl) > dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 & shuff_spk_t(:,kl) < end_time2(ik),kl);
    end    
    
    % convert to firing rate and store
    firing_rate(length(meta.TrialInfo_block1.trialT) + ik) = length(these_spk) / trial_block2.trial_length(ik);
    firing_rate_block2(ik) = length(these_spk) / trial_block2.trial_length(ik);
    for kl = 1:1000
    shuff_firing_rate(length(meta.TrialInfo_block1.trialT) + ik,kl) = length(these_spk_shuff{kl}) / trial_block2.trial_length(ik);
    end
    
end

end_time_trials = cat(1,end_time,end_time2);

FRATE.Task.Trial_firing_rate = firing_rate;
FRATE.Cue.Trial_firing_rate_block1 = firing_rate_block1;
FRATE.Cue.Trial_firing_rate_block2 = firing_rate_block2;
FRATE.Shuff.Trial_firing_rate_shuff = shuff_firing_rate;
FRATE.Interval.Trial = iv(start_time_trials,end_time_trials);
%% divided by rew vs unrew
rew_trial_num_block1 = 1;
unrew_trial_num_block1 = 1;
rew_trial_num_block2 = 1;
unrew_trial_num_block2 = 1;

for ik = 1:length(meta.TrialInfo_block1.trialT)
    switch meta.TrialInfo_block1.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_rew_trial(rew_trial_num_block1) = firing_rate(ik);
            fr_rew_block1(rew_trial_num_block1) = firing_rate(ik);
            rew_trial_num_block1 = rew_trial_num_block1 + 1;
        case 0
            fr_unrew_trial(unrew_trial_num_block1) = firing_rate(ik);
            fr_unrew_block1(unrew_trial_num_block1) = firing_rate(ik);
            unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
    end
end

for ik = 1:length(meta.TrialInfo_block2.trialT) - 1 %% -1 just to do trial minus last one as it isn't always completed
    switch meta.TrialInfo_block2.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_rew_trial(rew_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_rew_block2(rew_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            rew_trial_num_block1 = rew_trial_num_block1 + 1;
            rew_trial_num_block2 = rew_trial_num_block2 + 1;
        case 0
            fr_unrew_trial(unrew_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_unrew_block2(unrew_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
            unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
    end
end

FRATE.Reward.Trial_firing_rate_reward = fr_rew_trial;
FRATE.Reward.Trial_firing_rate_unreward = fr_unrew_trial;
FRATE.Condition.Trial_firing_rate_block1_rew = fr_rew_block1;
FRATE.Condition.Trial_firing_rate_block1_unrew = fr_unrew_block1;
FRATE.Condition.Trial_firing_rate_block2_rew = fr_rew_block2;
FRATE.Condition.Trial_firing_rate_block2_unrew = fr_unrew_block2;

%% divided by arm
arm1_trial_num_block1 = 1;
arm2_trial_num_block1 = 1;
arm3_trial_num_block1 = 1;
arm4_trial_num_block1 = 1;
arm1_trial_num_block2 = 1;
arm2_trial_num_block2 = 1;
arm3_trial_num_block2 = 1;
arm4_trial_num_block2 = 1;

for ik = 1:length(meta.TrialInfo_block1.trialT)
    switch meta.TrialInfo_block1.photosensorID(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_arm1_trial(arm1_trial_num_block1) = firing_rate(ik);
            fr_arm1_block1(arm1_trial_num_block1) = firing_rate(ik);
            arm1_trial_num_block1 = arm1_trial_num_block1 + 1;
        case 2
            fr_arm2_trial(arm2_trial_num_block1) = firing_rate(ik);
            fr_arm2_block1(arm2_trial_num_block1) = firing_rate(ik);
            arm2_trial_num_block1 = arm2_trial_num_block1 + 1;
        case 3
            fr_arm3_trial(arm3_trial_num_block1) = firing_rate(ik);
            fr_arm3_block1(arm3_trial_num_block1) = firing_rate(ik);
            arm3_trial_num_block1 = arm3_trial_num_block1 + 1;
        case 4
            fr_arm4_trial(arm4_trial_num_block1) = firing_rate(ik);
            fr_arm4_block1(arm4_trial_num_block1) = firing_rate(ik);
            arm4_trial_num_block1 = arm4_trial_num_block1 + 1;
    end
end

for ik = 1:length(meta.TrialInfo_block2.trialT)  - 1 %% -1 just to do trial minus last one as it isn't always completed
    switch meta.TrialInfo_block2.photosensorID(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_arm1_trial(arm1_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_arm1_block2(arm1_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            arm1_trial_num_block1 = arm1_trial_num_block1 + 1;
            arm1_trial_num_block2 = arm1_trial_num_block2 + 1;
        case 2
            fr_arm2_trial(arm2_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_arm2_block2(arm2_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            arm2_trial_num_block1 = arm2_trial_num_block1 + 1;
            arm2_trial_num_block2 = arm2_trial_num_block2 + 1;
        case 3
            fr_arm3_trial(arm3_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_arm3_block2(arm3_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            arm3_trial_num_block1 = arm3_trial_num_block1 + 1;
            arm3_trial_num_block2 = arm3_trial_num_block2 + 1;
        case 4
            fr_arm4_trial(arm4_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            fr_arm4_block2(arm4_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
            arm4_trial_num_block1 = arm4_trial_num_block1 + 1;
            arm4_trial_num_block2 = arm4_trial_num_block2 + 1;
    end
end

fr_arm1_trial = fr_arm1_trial';
fr_arm2_trial = fr_arm2_trial';
fr_arm3_trial = fr_arm3_trial';
fr_arm4_trial = fr_arm4_trial';
group1(1:length(fr_arm1_trial),1) = 1;
group2(1:length(fr_arm2_trial),1) = 2;
group3(1:length(fr_arm3_trial),1) = 3;
group4(1:length(fr_arm4_trial),1) = 4;
FRATE.Arm.Trial_firing_rate = cat(1,fr_arm1_trial,fr_arm2_trial,fr_arm3_trial,fr_arm4_trial);
FRATE.Arm.Trial_firing_rate_groups = cat(1,group1,group2,group3,group4);

%% app vs skip for unrew
skip_trial_num_block1 = 1;
app_trial_num_block1 = 1;
skip_trial_num_block2 = 1;
app_trial_num_block2 = 1;

for ik = 1:length(meta.TrialInfo_block1.trialT)
    if meta.TrialInfo_block1.rewarded(ik) == 0
        switch meta.TrialInfo_block1.nosepokeID(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 0
                fr_skip(skip_trial_num_block1) = firing_rate(ik);
                fr_skip_block1(skip_trial_num_block1) = firing_rate(ik);
                skip_trial_num_block1 = skip_trial_num_block1 + 1;
            case {5,6,7,8}
                fr_app(app_trial_num_block1) = firing_rate(ik);
                fr_app_block1(app_trial_num_block1) = firing_rate(ik);
                app_trial_num_block1 = app_trial_num_block1 + 1;
        end
    end
end

for ik = 1:length(meta.TrialInfo_block2.trialT)  - 1 %% -1 just to do trial minus last one as it isn't always completed
    if meta.TrialInfo_block2.rewarded(ik) == 0
        switch meta.TrialInfo_block2.nosepokeID(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 0
                fr_skip(skip_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
                fr_skip_block2(skip_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
                skip_trial_num_block1 = skip_trial_num_block1 + 1;
                skip_trial_num_block2 = skip_trial_num_block2 + 1;
            case {5,6,7,8}
                fr_app(app_trial_num_block1) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
                fr_app_block2(app_trial_num_block2) = firing_rate(length(meta.TrialInfo_block1.trialT) + ik);
                app_trial_num_block1 = app_trial_num_block1 + 1;
                app_trial_num_block2 = app_trial_num_block2 + 1;
        end
    end
end

FRATE.Approach.Trial_firing_rate_skip = fr_skip;
FRATE.Approach.Trial_firing_rate_app = fr_app;

%% 1 s before trial

end_time = nan(length(meta.TrialInfo_block1.trialT),1);
end_time2 = nan(length(meta.TrialInfo_block2.trialT),1);
start_time = [];
start_time2 = [];
firing_rate_b4_trial = cat(1,end_time,end_time2);
firing_rate_b4_trial_block1 = nan(length(meta.TrialInfo_block1.trialT),1);
firing_rate_b4_trial_block2 = nan(length(meta.TrialInfo_block2.trialT),1);

for ik = 1:length(meta.TrialInfo_block1.trialT)
    end_time(ik) = dataPoint.Trials(ik)/1000;

    if meta.TrialInfo_block1.trial_length_analysis(ik) < 1
        start_time(ik) = dataPoint.Trials(ik)/1000 - meta.TrialInfo_block1.trial_length_analysis(ik);
    else
        start_time(ik) = dataPoint.Trials(ik)/1000 - 1;
    end
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > start_time(ik) & spk_t < end_time(ik));
    for kl = 1:1000
        these_spk_shuff{kl} = shuff_spk_t(shuff_spk_t(:,kl) > start_time(ik) & shuff_spk_t(:,kl) < end_time(ik),kl);
    end
    
    % convert to firing rate and store
    firing_rate_b4_trial(ik) = length(these_spk) / 1;
    firing_rate_b4_trial_block1(ik) = length(these_spk) / 1;
    for kl = 1:1000
        shuff_firing_rate_b4_trial(ik,kl) = length(these_spk_shuff{kl}) / 1;
    end
end

for ik = 1:length(meta.TrialInfo_block2.trialT) %-1 % -1 just to do trial minus last one as it isn't always completed
    end_time2(ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000;
    
    if meta.TrialInfo_block2.trial_length_analysis(ik) < 1
        start_time2(ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 - meta.TrialInfo_block2.trial_length_analysis(ik);
    else
        start_time2(ik) = dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000 - 1;
    end
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > start_time2(ik) & spk_t < end_time2(ik));
    for kl = 1:1000
        these_spk_shuff{kl} = shuff_spk_t(shuff_spk_t(:,kl) > start_time2(ik) & shuff_spk_t(:,kl) < end_time2(ik),kl);
    end
    
    % convert to firing rate and store
    firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik) = length(these_spk) / 1;
    firing_rate_b4_trial_block2(ik) = length(these_spk) / 1;
    for kl = 1:1000
        shuff_firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik,kl) = length(these_spk_shuff{kl}) / 1;
    end
end

FRATE.Task.Trial_B4_firing_rate = firing_rate_b4_trial;
FRATE.Cue.Trial_B4_firing_rate_block1 = firing_rate_b4_trial_block1;
FRATE.Cue.Trial_B4_firing_rate_block2 = firing_rate_b4_trial_block2;
FRATE.Shuff.Trial_B4_firing_rate_shuff = shuff_firing_rate_b4_trial;

%% rew vs unrew
rew_trial_num_block1 = 1;
unrew_trial_num_block1 = 1;
rew_trial_num_block2 = 1;
unrew_trial_num_block2 = 1;

for ik = 2:length(meta.TrialInfo_block1.trialT)
    switch meta.TrialInfo_block1.rewarded(ik - 1) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_B4_trial_rew_block1(rew_trial_num_block1) = firing_rate_b4_trial(ik);
            rew_trial_num_block1 = rew_trial_num_block1 + 1;
        case 0
            fr_B4_trial_unrew_block1(unrew_trial_num_block1) = firing_rate_b4_trial(ik);
            unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
    end
end

for ik = 2:length(meta.TrialInfo_block2.trialT)  - 1 %% -1 just to do trial minus last one as it isn't always completed
    switch meta.TrialInfo_block2.rewarded(ik - 1) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_B4_trial_rew_block2(rew_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            rew_trial_num_block2 = rew_trial_num_block2 + 1;
        case 0
            fr_B4_trial_unrew_block2(unrew_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
    end
end

FRATE.Condition.Trial_B4_firing_rate_block1_rew = fr_B4_trial_rew_block1;
FRATE.Condition.Trial_B4_firing_rate_block1_unrew =  fr_B4_trial_unrew_block1;
FRATE.Condition.Trial_B4_firing_rate_block2_rew = fr_B4_trial_rew_block2;
FRATE.Condition.Trial_B4_firing_rate_block2_unrew =  fr_B4_trial_unrew_block2;
%% divided by arm
arm1_trial_num_block1 = 1;
arm2_trial_num_block1 = 1;
arm3_trial_num_block1 = 1;
arm4_trial_num_block1 = 1;
arm1_trial_num_block2 = 1;
arm2_trial_num_block2 = 1;
arm3_trial_num_block2 = 1;
arm4_trial_num_block2 = 1;

for ik = 2:length(meta.TrialInfo_block1.trialT)
    switch meta.TrialInfo_block1.photosensorID(ik - 1) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_B4_trial_arm1(arm1_trial_num_block1) = firing_rate_b4_trial(ik);
            fr_B4_trial_arm1_block1(arm1_trial_num_block1) = firing_rate_b4_trial(ik);
            arm1_trial_num_block1 = arm1_trial_num_block1 + 1;
        case 2
            fr_B4_trial_arm2(arm2_trial_num_block1) = firing_rate_b4_trial(ik);
            fr_B4_trial_arm2_block1(arm2_trial_num_block1) = firing_rate_b4_trial(ik);
            arm2_trial_num_block1 = arm2_trial_num_block1 + 1;
        case 3
            fr_B4_trial_arm3(arm3_trial_num_block1) = firing_rate_b4_trial(ik);
            fr_B4_trial_arm3_block1(arm3_trial_num_block1) = firing_rate_b4_trial(ik);
            arm3_trial_num_block1 = arm3_trial_num_block1 + 1;
        case 4
            fr_B4_trial_arm4(arm4_trial_num_block1) = firing_rate_b4_trial(ik);
            fr_B4_trial_arm4_block1(arm4_trial_num_block1) = firing_rate_b4_trial(ik);
            arm4_trial_num_block1 = arm4_trial_num_block1 + 1;
    end
end

for ik = 2:length(meta.TrialInfo_block2.trialT)  - 1 %% -1 just to do trial minus last one as it isn't always completed
    switch meta.TrialInfo_block2.photosensorID(ik - 1) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            fr_B4_trial_arm1(arm1_trial_num_block1) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            fr_B4_trial_arm1_block2(arm1_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            arm1_trial_num_block1 = arm1_trial_num_block1 + 1;
            arm1_trial_num_block2 = arm1_trial_num_block2 + 1;
        case 2
            fr_B4_trial_arm2(arm2_trial_num_block1) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            fr_B4_trial_arm2_block2(arm2_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            arm2_trial_num_block1 = arm2_trial_num_block1 + 1;
            arm2_trial_num_block2 = arm2_trial_num_block2 + 1;
        case 3
            fr_B4_trial_arm3(arm3_trial_num_block1) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            fr_B4_trial_arm3_block2(arm3_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            arm3_trial_num_block1 = arm3_trial_num_block1 + 1;
            arm3_trial_num_block2 = arm3_trial_num_block2 + 1;
        case 4
            fr_B4_trial_arm4(arm4_trial_num_block1) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            fr_B4_trial_arm4_block2(arm4_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
            arm4_trial_num_block1 = arm4_trial_num_block1 + 1;
            arm4_trial_num_block2 = arm4_trial_num_block2 + 1;
    end
end

fr_B4_trial_arm1 = fr_B4_trial_arm1';
fr_B4_trial_arm2 = fr_B4_trial_arm2';
fr_B4_trial_arm3 = fr_B4_trial_arm3';
fr_B4_trial_arm4 = fr_B4_trial_arm4';
group_B4_1(1:length(fr_B4_trial_arm1),1) = 1;
group_B4_2(1:length(fr_B4_trial_arm2),1) = 2;
group_B4_3(1:length(fr_B4_trial_arm3),1) = 3;
group_B4_4(1:length(fr_B4_trial_arm4),1) = 4;
FRATE.Arm.Trial_B4_firing_rate = cat(1,fr_B4_trial_arm1,fr_B4_trial_arm2,fr_B4_trial_arm3,fr_B4_trial_arm4);
FRATE.Arm.Trial_B4_firing_rate_groups = cat(1,group_B4_1,group_B4_2,group_B4_3,group_B4_4);

%% unrew app vs skip
skip_trial_num_block1 = 1;
app_trial_num_block1 = 1;
skip_trial_num_block2 = 1;
app_trial_num_block2 = 1;

for ik = 2:length(meta.TrialInfo_block1.trialT)
    if meta.TrialInfo_block1.rewarded(ik - 1) == 0
        switch meta.TrialInfo_block1.nosepokeID(ik - 1) %add to skiparded or apparded count depending on if trial was skiparded or not
            case 0
                fr_B4_trial_skip_block1(skip_trial_num_block1) = firing_rate_b4_trial(ik);
                skip_trial_num_block1 = skip_trial_num_block1 + 1;
            case {5,6,7,8}
                fr_B4_trial_app_block1(app_trial_num_block1) = firing_rate_b4_trial(ik);
                app_trial_num_block1 = app_trial_num_block1 + 1;
        end
    end
end

for ik = 2:length(meta.TrialInfo_block2.trialT)  - 1 %% -1 just to do trial minus last one as it isn't always completed
    if meta.TrialInfo_block2.rewarded(ik - 1) == 0
        switch meta.TrialInfo_block2.nosepokeID(ik - 1) %add to rewarded or apparded count depending on if trial was rewarded or not
            case 0
                fr_B4_trial_skip_block2(skip_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
                skip_trial_num_block2 = skip_trial_num_block2 + 1;
            case {5,6,7,8}
                fr_B4_trial_app_block2(app_trial_num_block2) = firing_rate_b4_trial(length(meta.TrialInfo_block1.trialT) + ik);
                app_trial_num_block2 = app_trial_num_block2 + 1;
        end
    end
end

%% data TRIALS during NOSEPOKES
nosepoke_number_trial = 1;
nosepoke_number_trial2 = 1;
end_time = nan(length(meta.TrialInfo_block1.trialT),1);
end_time2 = nan(length(meta.TrialInfo_block2.trialT),1);
start_time = [];
start_time2 = [];
% firing_rate_trial_np = cat(1,end_time,end_time2);
% firing_rate_trial_np_block1 = nan(length(meta.TrialInfo_block1.trialT),1);
% firing_rate_trial_np_block2 = nan(length(meta.TrialInfo_block2.trialT),1);

for ik = 1:length(meta.TrialInfo_block1.trialT)
    if meta.TrialInfo_block1.nosepokeT(ik) ~= 0
        if meta.TrialInfo_block1.nosepoke_length(ik) < 1
            trial_np_block1.trial_length(nosepoke_number_trial) = meta.TrialInfo_block1.nosepoke_length(ik);
        else
            trial_np_block1.trial_length(nosepoke_number_trial) = 1;
        end
        end_time(nosepoke_number_trial) = dataPoint.Nosepokes(nosepoke_number_trial)/1000;
        start_time(nosepoke_number_trial) = dataPoint.Nosepokes(nosepoke_number_trial)/1000 - trial_np_block1.trial_length(nosepoke_number_trial);
        
        % now, count spikes between start and end
        these_spk = spk_t(spk_t > start_time(nosepoke_number_trial) & spk_t < end_time(nosepoke_number_trial));
        for kl = 1:1000
        these_spk_shuff{kl} = shuff_spk_t(shuff_spk_t(:,kl) > start_time(nosepoke_number_trial) & shuff_spk_t(:,kl) < end_time(nosepoke_number_trial),kl);
        end
        
        % convert to firing rate and store
        firing_rate_trial_np(nosepoke_number_trial) = length(these_spk) / trial_np_block1.trial_length(nosepoke_number_trial);
        firing_rate_trial_np_block1(nosepoke_number_trial) = length(these_spk) / trial_np_block1.trial_length(nosepoke_number_trial);
        for kl = 1:1000
        shuff_firing_rate_trial_np(nosepoke_number_trial,kl) = length(these_spk_shuff{kl}) / trial_np_block1.trial_length(nosepoke_number_trial);
    end
        nosepoke_number_trial = nosepoke_number_trial + 1;
    end
end

for ik = 1:length(meta.TrialInfo_block2.trialT) %-1 % -1 just to do trial minus last one as it isn't always completed
    if meta.TrialInfo_block2.nosepokeT(ik) ~= 0
        if meta.TrialInfo_block2.nosepoke_length(ik) < 1
            trial_np_block2.trial_length(nosepoke_number_trial2) = meta.TrialInfo_block2.nosepoke_length(ik);
        else
            trial_np_block2.trial_length(nosepoke_number_trial2) = 1;
        end
        end_time2(nosepoke_number_trial2) = dataPoint.Nosepokes(nosepoke_number_trial)/1000;
        start_time2(nosepoke_number_trial2) = dataPoint.Nosepokes(nosepoke_number_trial)/1000 - trial_np_block2.trial_length(nosepoke_number_trial2);
        
        % now, count spikes between start and end
        these_spk = spk_t(spk_t > start_time2(nosepoke_number_trial2) & spk_t < end_time2(nosepoke_number_trial2));
        for kl = 1:1000
        these_spk_shuff{kl} = shuff_spk_t(shuff_spk_t(:,kl) > start_time2(nosepoke_number_trial2) & shuff_spk_t(:,kl) < end_time2(nosepoke_number_trial2),kl);
        end
        
        % convert to firing rate and store
        firing_rate_trial_np(nosepoke_number_trial) = length(these_spk) / trial_np_block2.trial_length(nosepoke_number_trial2);
        firing_rate_trial_np_block2(nosepoke_number_trial2) = length(these_spk) / trial_np_block2.trial_length(nosepoke_number_trial2);
        for kl = 1:1000
        shuff_firing_rate_trial_np(nosepoke_number_trial,kl) = length(these_spk_shuff{kl}) / trial_np_block2.trial_length(nosepoke_number_trial2);
    end
        nosepoke_number_trial = nosepoke_number_trial + 1;
        nosepoke_number_trial2 = nosepoke_number_trial2 + 1;
    end
end

FRATE.Task.Trial_np_firing_rate = firing_rate_trial_np;
FRATE.Shuff.Trial_np_firing_rate = shuff_firing_rate_trial_np;
toc