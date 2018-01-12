function RAST = genRAST(meta,spk_t,dataPoint)
% function RAST = genRAST(meta,spk_t,dataPoint)
%
% computes average per-event time histogram of spike counts relative to
% specified events
%
% INPUTS:
% meta: meta file containing AMPX and task data
% spk_t: vector of spike times
% dataPoint: list of trial start times (in ~ms)
%
% OUTPUTS:
% RAST: struct of spikes by trial for rasterplot

%% for trials
disp('generating rasterplots')

for ik = 1:length(meta.TrialInfo_block1.trialT)
    if isempty (find(spk_t > dataPoint.Trials(ik)/1000-3 & spk_t < dataPoint.Trials(ik)/1000+6));
        trial_spikes_block1{ik} = 0;
        trial_spikes_time{ik} = 0;
    else
        trial_spikes_block1{ik} = find(spk_t > dataPoint.Trials(ik)/1000-3 & spk_t < dataPoint.Trials(ik)/1000+6);
        trial_spikes_time{ik} = spk_t(trial_spikes_block1{ik}) - dataPoint.Trials(ik)/1000;
    end  
end

for ik = 1:length(meta.TrialInfo_block2.trialT)
    if isempty (find(spk_t > dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000-3 & spk_t < dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000+6));
        trial_spikes_block2{ik} = 0;
        trial_spikes_time{length(meta.TrialInfo_block1.trialT) + ik} = 0;
    else
        trial_spikes_block2{ik} = find(spk_t > dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000-3 & spk_t < dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000+6);
        trial_spikes_time{length(meta.TrialInfo_block1.trialT) + ik} = spk_t(trial_spikes_block2{ik}) - dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ik)/1000;
    end  
end

RAST.Trial = trial_spikes_time;