%% synthetic spike train modulated by cue-onset for reward-availablility
TrialInfo_b1 = {'TrialInfo_block1' 'TrialInfo{1,1}'};
TrialInfo_b2 = {'TrialInfo_block2' 'TrialInfo{1,2}'};
AnalType = {'rewarded_cue' 'cue_location'};
Trial_L = {'trial_length_analysis' 'trial_length'};
% Trial_type = 1;

mat_files = dir('*.mat');
for iCell = 352:length(dir('*.mat'))
    disp(iCell);
    load(mat_files(iCell).name);
    switch strcmp(mat_files(iCell).name(1:4),'R060');
    case 0
        old_v_new = 1;
        case 1
            old_v_new = 2;
    end
    switch mod(iCell,2)
        case 0
            AType = 1;
            Trial_type = 2;
        case 1
            AType = 2;
            Trial_type = 1;
    end
    synth.rewarded_cue{iCell} = cat(1,metadata.TrialInfo{1,1}.rewarded,metadata.TrialInfo{1,2}.rewarded);
    synth.cue_location{iCell} = cat(1,metadata.TrialInfo{1,1}.photosensorID,metadata.TrialInfo{1,2}.photosensorID);
    synth.trial_length{iCell} = cat(1,metadata.TrialInfo{1,1}.(Trial_L{Trial_type}),metadata.TrialInfo{1,2}.(Trial_L{Trial_type}));
    synth.np_length{iCell} = cat(1,metadata.TrialInfo{1,1}.nosepoke_length,metadata.TrialInfo{1,2}.nosepoke_length);
    spk_count = 1;
    for iTrial = 1:length(dataPoint.Trials)
        if synth.(AnalType{AType}){iCell}(iTrial) == 1;
        iTime = .05;
        switch AType
            case 1
        t_length = synth.trial_length{iCell}(iTrial) + synth.np_length{iCell}(iTrial);
            case 2
                t_length = synth.trial_length{iCell}(iTrial);
        end
                
        while iTime < t_length
   synth.spk_t{iCell}(spk_count) = dataPoint.Trials(iTrial)/1000 + iTime;
   iTime = iTime + .05;
   spk_count = spk_count + 1;
        end
        end
    end
end