function dataPoint = genDP(meta)
% function dataPoint = genDP(meta)
%
% Finds the times of cue-onset for the neural data from the MATLAB record
% of the behavior
%
% INPUTS:
% meta: ts of cue-onsets for behavioral data
%
% OUTPUTS:
% dataPoint: ts of cue-onsets for neural data

%%
tic
disp('getting event times')
t = [0 4000];
binsize = 0.001;
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
%% Generate datapoints for trials
for i = 1:length(meta.TrialInfo_block1.trialT)
    ZeroT_trials(i) = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.TrialInfo_block1.trialT(i) + meta.offsetT_block1(i); %apply offset to MATLAB data
    dataPoint.Trials(i) = find(tbin_centers > ZeroT_trials(i), 1, 'first'); %find times for all events
end % block 1
for j = 1:length(meta.TrialInfo_block2.trialT)
    ZeroT_trials(length(meta.TrialInfo_block1.trialT)  + j) = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.TrialInfo_block2.trialT(j) + meta.offsetT_block2(j); %apply offset to MATLAB data
    dataPoint.Trials(length(meta.TrialInfo_block1.trialT)  + j) = find(tbin_centers > ZeroT_trials(length(meta.TrialInfo_block1.trialT)  + j), 1, 'first'); %find times for all events
end % add block 2 to block 1
%% Generate datapoints for nosepokes
nosepoke_counter = 1;
for i = 1:length(meta.TrialInfo_block1.nosepokeT)
    if meta.TrialInfo_block1.nosepokeT(i) ~= 0
        ZeroT_nosepokes(i) = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.TrialInfo_block1.nosepokeT(i) + meta.offsetT_block1(i); %apply offset to MATLAB data
        dataPoint.Nosepokes(nosepoke_counter) = find(tbin_centers > ZeroT_nosepokes(i), 1, 'first'); %find times for all events
        nosepoke_counter = nosepoke_counter + 1;
    end
end
block2_start = nosepoke_counter - 1;
for j = 1:length(meta.TrialInfo_block2.nosepokeT)
    if meta.TrialInfo_block2.nosepokeT(j) ~= 0
        ZeroT_nosepokes(block2_start  + j) = meta.data_pre.tvec(end) + meta.gap1.tvec(end) + meta.data_block1.tvec(end) + meta.gap2.tvec(end) + meta.TrialInfo_block2.nosepokeT(j) + meta.offsetT_block2(j); %apply offset to MATLAB data
        dataPoint.Nosepokes(nosepoke_counter) = find(tbin_centers > ZeroT_nosepokes(block2_start  + j), 1, 'first'); %find times for all events
        nosepoke_counter = nosepoke_counter + 1;
    end
end % add block 2 to block 1
toc