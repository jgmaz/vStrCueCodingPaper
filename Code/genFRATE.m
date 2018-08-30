function FRATE = genFRATE(meta,spk_t)
% function FRATE = genFRATE(meta,spk_t)
%
% Generates average, and trial-by-trial firing rates broken down a variety
% of ways
%
% INPUTS:
% meta: meta file containing AMPX and task data
% spk_t: vector of spike times
% meta.dataPoint: list of trial and nosepoke start times (in ~ms)
%
% OUTPUTS:
% FRATE: struct containing pre- and post-trial firing rates

%% trials
end_time = nan(length(meta.TrialInfo{1,1}.trialT),1);
end_time2 = nan(length(meta.TrialInfo{1,2}.trialT),1);
firing_rate = cat(1,end_time,end_time2);
firing_rate_block1 = nan(length(meta.TrialInfo{1,1}.trialT),1);
firing_rate_block2 = nan(length(meta.TrialInfo{1,2}.trialT),1);

for ik = 1:length(meta.TrialInfo{1,1}.trialT)
    if meta.TrialInfo{1,1}.trial_length_analysis(ik) < 1
        trial_block1.trial_length(ik) = meta.TrialInfo{1,1}.trial_length_analysis(ik);
        start_time_trials(ik) = meta.dataPoint.Trials(ik)/1000 - meta.TrialInfo{1,1}.trial_length_analysis(ik);
    else
        trial_block1.trial_length(ik) = 1;
        start_time_trials(ik) = meta.dataPoint.Trials(ik)/1000 - 1;
    end
    end_time(ik) = meta.dataPoint.Trials(ik)/1000 + trial_block1.trial_length(ik);
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > meta.dataPoint.Trials(ik)/1000 & spk_t < end_time(ik));
    
    % convert to firing rate and store
    firing_rate(ik) = length(these_spk) / trial_block1.trial_length(ik);
    firing_rate_block1(ik) = length(these_spk) / trial_block1.trial_length(ik);    
end

for ik = 1:length(meta.TrialInfo{1,2}.trialT) %-1 % -1 just to do trial minus last one as it isn't always completed
    if meta.TrialInfo{1,2}.trial_length_analysis(ik) < 1
        trial_block2.trial_length(ik) = meta.TrialInfo{1,2}.trial_length_analysis(ik);
        start_time_trials(length(meta.TrialInfo{1,1}.trialT) + ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 - meta.TrialInfo{1,2}.trial_length_analysis(ik);
    else
        trial_block2.trial_length(ik) = 1;
        start_time_trials(length(meta.TrialInfo{1,1}.trialT) + ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 - 1;
    end
    end_time2(ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 + trial_block2.trial_length(ik);
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 & spk_t < end_time2(ik));
    % convert to firing rate and store
    firing_rate(length(meta.TrialInfo{1,1}.trialT) + ik) = length(these_spk) / trial_block2.trial_length(ik);
    firing_rate_block2(ik) = length(these_spk) / trial_block2.trial_length(ik);
    
end

end_time_trials = cat(1,end_time,end_time2);

FRATE.Task.Trial_firing_rate = firing_rate;
FRATE.Cue.Trial_firing_rate_block1 = firing_rate_block1;
FRATE.Cue.Trial_firing_rate_block2 = firing_rate_block2;
FRATE.Interval.Trial = iv(start_time_trials,end_time_trials);

%% 1 s before trial

end_time = nan(length(meta.TrialInfo{1,1}.trialT),1);
end_time2 = nan(length(meta.TrialInfo{1,2}.trialT),1);
start_time = [];
start_time2 = [];
firing_rate_b4_trial = cat(1,end_time,end_time2);
firing_rate_b4_trial_block1 = nan(length(meta.TrialInfo{1,1}.trialT),1);
firing_rate_b4_trial_block2 = nan(length(meta.TrialInfo{1,2}.trialT),1);

for ik = 1:length(meta.TrialInfo{1,1}.trialT)
    end_time(ik) = meta.dataPoint.Trials(ik)/1000;

    if meta.TrialInfo{1,1}.trial_length_analysis(ik) < 1
        start_time(ik) = meta.dataPoint.Trials(ik)/1000 - meta.TrialInfo{1,1}.trial_length_analysis(ik);
    else
        start_time(ik) = meta.dataPoint.Trials(ik)/1000 - 1;
    end
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > start_time(ik) & spk_t < end_time(ik));
    
    % convert to firing rate and store
    firing_rate_b4_trial(ik) = length(these_spk) / 1;
    firing_rate_b4_trial_block1(ik) = length(these_spk) / 1;
end

for ik = 1:length(meta.TrialInfo{1,2}.trialT) %-1 % -1 just to do trial minus last one as it isn't always completed
    end_time2(ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000;
    
    if meta.TrialInfo{1,2}.trial_length_analysis(ik) < 1
        start_time2(ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 - meta.TrialInfo{1,2}.trial_length_analysis(ik);
    else
        start_time2(ik) = meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000 - 1;
    end
    
    % now, count spikes between start and end
    these_spk = spk_t(spk_t > start_time2(ik) & spk_t < end_time2(ik));

    % convert to firing rate and store
    firing_rate_b4_trial(length(meta.TrialInfo{1,1}.trialT) + ik) = length(these_spk) / 1;
    firing_rate_b4_trial_block2(ik) = length(these_spk) / 1;

end

FRATE.Task.Trial_B4_firing_rate = firing_rate_b4_trial;
FRATE.Cue.Trial_B4_firing_rate_block1 = firing_rate_b4_trial_block1;
FRATE.Cue.Trial_B4_firing_rate_block2 = firing_rate_b4_trial_block2;