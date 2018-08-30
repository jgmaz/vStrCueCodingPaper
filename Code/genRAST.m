function RAST = genRAST(meta,spk_t)
% function RAST = genRAST(meta,spk_t)
%
% computes average per-event time histogram of spike counts relative to
% specified events
%
% INPUTS:
% meta: meta file containing AMPX and task data
% spk_t: vector of spike times
% meta.dataPoint: list of trial start times (in ~ms)
%
% OUTPUTS:
% RAST: struct of spikes by trial for rasterplot

%% for trials
disp('generating rasterplots')

for ik = 1:length(meta.TrialInfo{1,1}.trialT)
    if isempty (find(spk_t > meta.dataPoint.Trials(ik)/1000-3 & spk_t < meta.dataPoint.Trials(ik)/1000+6));
        trial_spikes_block1{ik} = 0;
        trial_spikes_time{ik} = 0;
    else
        trial_spikes_block1{ik} = find(spk_t > meta.dataPoint.Trials(ik)/1000-3 & spk_t < meta.dataPoint.Trials(ik)/1000+6);
        trial_spikes_time{ik} = spk_t(trial_spikes_block1{ik}) - meta.dataPoint.Trials(ik)/1000;
    end  
end

for ik = 1:length(meta.TrialInfo{1,2}.trialT)
    if isempty (find(spk_t > meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000-3 & spk_t < meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000+6));
        trial_spikes_block2{ik} = 0;
        trial_spikes_time{length(meta.TrialInfo{1,1}.trialT) + ik} = 0;
    else
        trial_spikes_block2{ik} = find(spk_t > meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000-3 & spk_t < meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000+6);
        trial_spikes_time{length(meta.TrialInfo{1,1}.trialT) + ik} = spk_t(trial_spikes_block2{ik}) - meta.dataPoint.Trials(length(meta.TrialInfo{1,1}.trialT) + ik)/1000;
    end  
end

RAST.Trial = trial_spikes_time;