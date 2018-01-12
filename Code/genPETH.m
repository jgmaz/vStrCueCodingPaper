function PETH = genPETH(sesh,meta,spk_t,dataPoint)
% function PETH = genPETH(sesh,meta,spk_t,dataPoint)
%
% computes average peri-event time histogram of spike counts relative to
% specified events
%
% INPUTS:
% sesh: variable containing session specific information inputed in
% overview_graph
% meta: meta file containing AMPX and task data
% spk_t: vector of spike times
% dataPoint: list of trial start times (in ~ms)
%
% OUTPUTS:
% PETH.Trial: peri-event time histograms according to trial start by cue
% PETH.Arm: peri-event time histograms according to trial start by location
% PETH.Approach: peri-event time histograms according to trial start for unrewarded cues separated by whether or not they
% approached the receptacle (error trials)
% PETH.Trial_np: peri-event time histograms according to trial start by cue
% for only approach trials
%
% CONFIGS:
% sesh.PETH = generate a PETH for each event type specfied with a '1'

%% SDF
tic
disp('generating peri-event histograms')

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

%% rewarded vs unrewarded trials for both blocks
if sesh.PETH.Trial == 1
    
    temp = -5000;
    rew_trial_num_block1 = 1;
    unrew_trial_num_block1 = 1;
    rew_trial_num_block2 = 1;
    unrew_trial_num_block2 = 1;
    trial_num = 1;
    trial_num2 = 1;
    for ik = 1:length(meta.TrialInfo_block1.trialT)
        switch meta.TrialInfo_block1.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                for jk = 1:15001
                    if ik == 1
                        rew_trials_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    rew_trials_block1(jk,rew_trial_num_block1) = NaN;
                                    trials_PETH(jk,trial_num) = NaN;
                                    trials_block1_PETH(jk,trial_num) = NaN;
                                else
                                    rew_trials_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                rew_trials_block1(jk,rew_trial_num_block1) = NaN;
                                trials_PETH(jk,trial_num) = NaN;
                                trials_block1_PETH(jk,trial_num) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            rew_trials_block1(jk:15001,rew_trial_num_block1) = NaN;
                            trials_PETH(jk:15001,trial_num) = NaN;
                            trials_block1_PETH(jk:15001,trial_num) = NaN;
                            break
                        else
                            rew_trials_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                rew_trial_num_block1 = rew_trial_num_block1 + 1;
                trial_num = trial_num + 1;
            case 0
                for jk = 1:15001
                    if ik == 1
                        unrew_trials_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    unrew_trials_block1(jk,unrew_trial_num_block1) = NaN;
                                    trials_PETH(jk,trial_num) = NaN;
                                    trials_block1_PETH(jk,trial_num) = NaN;
                                else
                                    unrew_trials_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                unrew_trials_block1(jk,unrew_trial_num_block1) = NaN;
                                trials_PETH(jk,trial_num) = NaN;
                                trials_block1_PETH(jk,trial_num) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            unrew_trials_block1(jk:15001,unrew_trial_num_block1) = NaN;
                            trials_PETH(jk:15001,trial_num) = NaN;
                            trials_block1_PETH(jk:15001,trial_num) = NaN;
                            break
                        else
                            unrew_trials_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
                trial_num = trial_num + 1;
        end
        
        temp = -5000;
    end
    
    for ip = 1:length(meta.TrialInfo_block2.trialT)
        switch meta.TrialInfo_block2.rewarded(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                for jk = 1:15001
                    if ip == 1
                        rew_trials_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    rew_trials_block2(jk,rew_trial_num_block2) = NaN;
                                    trials_PETH(jk,trial_num) = NaN;
                                    trials_block2_PETH(jk,trial_num2) = NaN;
                                else
                                    rew_trials_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                rew_trials_block2(jk,rew_trial_num_block2) = NaN;
                                trials_PETH(jk,trial_num) = NaN;
                                trials_block2_PETH(jk,trial_num2) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            rew_trials_block2(jk:15001,rew_trial_num_block2) = NaN;
                            trials_PETH(jk:15001,trial_num) = NaN;
                            trials_block2_PETH(jk:15001,trial_num2) = NaN;
                            break
                        else
                            rew_trials_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                rew_trial_num_block2 = rew_trial_num_block2 + 1;
                trial_num = trial_num + 1;
                trial_num2 = trial_num2 + 1;
            case 0
                for jk = 1:15001
                    if ip == 1
                        unrew_trials_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    unrew_trials_block2(jk,unrew_trial_num_block2) = NaN;
                                    trials_PETH(jk,trial_num) = NaN;
                                    trials_block2_PETH(jk,trial_num2) = NaN;
                                else
                                    unrew_trials_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                unrew_trials_block2(jk,unrew_trial_num_block2) = NaN;
                                trials_PETH(jk,trial_num) = NaN;
                                trials_block2_PETH(jk,trial_num2) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            unrew_trials_block2(jk:15001,unrew_trial_num_block2) = NaN;
                            trials_PETH(jk:15001,trial_num) = NaN;
                            trials_block2_PETH(jk:15001,trial_num2) = NaN;
                            break
                        else
                            unrew_trials_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
                trial_num = trial_num + 1;
                trial_num2 = trial_num2 + 1;
        end
        temp = -5000;
    end
    
    trials_rew_PETH = cat(2,rew_trials_block1,rew_trials_block2);
    trials_unrew_PETH = cat(2,unrew_trials_block1,unrew_trials_block2);
    
    for jj = 1:15001 %generate averaged responses
        if jj <= length(rew_trials_block1)
            avg_rew_trials_block1(jj,1) = nanmean(rew_trials_block1(jj,:));
            SEM_rew_trials_block1(jj,1) = nanstd(rew_trials_block1(jj,:)/sqrt(numel(rew_trials_block1(jj,:))-sum(isnan(rew_trials_block1(jj,:)))));
        end
        if jj <= length(unrew_trials_block1)
            avg_unrew_trials_block1(jj,1) = nanmean(unrew_trials_block1(jj,:));
            SEM_unrew_trials_block1(jj,1) = nanstd(unrew_trials_block1(jj,:)/sqrt(numel(unrew_trials_block1(jj,:))-sum(isnan(unrew_trials_block1(jj,:)))));
        end
        if jj <= length(rew_trials_block2)
            avg_rew_trials_block2(jj,1) = nanmean(rew_trials_block2(jj,:));
            SEM_rew_trials_block2(jj,1) = nanstd(rew_trials_block2(jj,:)/sqrt(numel(rew_trials_block2(jj,:))-sum(isnan(rew_trials_block2(jj,:)))));
        end
        if jj <= length(unrew_trials_block2)
            avg_unrew_trials_block2(jj,1) = nanmean(unrew_trials_block2(jj,:));
            SEM_unrew_trials_block2(jj,1) = nanstd(unrew_trials_block2(jj,:)/sqrt(numel(unrew_trials_block2(jj,:))-sum(isnan(unrew_trials_block2(jj,:)))));
        end
        
        if jj <= length(trials_rew_PETH)
            avg_trials_rew_PETH(jj,1) = nanmean(trials_rew_PETH(jj,:));
            SEM_trials_rew_PETH(jj,1) = nanstd(trials_rew_PETH(jj,:)/sqrt(numel(trials_rew_PETH(jj,:))-sum(isnan(trials_rew_PETH(jj,:)))));
        end
        
        if jj <= length(trials_unrew_PETH)
            avg_trials_unrew_PETH(jj,1) = nanmean(trials_unrew_PETH(jj,:));
            SEM_trials_unrew_PETH(jj,1) = nanstd(trials_unrew_PETH(jj,:)/sqrt(numel(trials_unrew_PETH(jj,:))-sum(isnan(trials_unrew_PETH(jj,:)))));
        end
        
        if jj <= length(trials_PETH)
            avg_trials_PETH(jj,1) = nanmean(trials_PETH(jj,:));
            SEM_trials_PETH(jj,1) = nanstd(trials_PETH(jj,:)/sqrt(numel(trials_PETH(jj,:))-sum(isnan(trials_PETH(jj,:)))));
        end
        
        if jj <= length(trials_block1_PETH)
            avg_trials_block1_PETH(jj,1) = nanmean(trials_block1_PETH(jj,:));
            SEM_trials_block1_PETH(jj,1) = nanstd(trials_block1_PETH(jj,:)/sqrt(numel(trials_block1_PETH(jj,:))-sum(isnan(trials_block1_PETH(jj,:)))));
        end
        
        if jj <= length(trials_block2_PETH)
            avg_trials_block2_PETH(jj,1) = nanmean(trials_block2_PETH(jj,:));
            SEM_trials_block2_PETH(jj,1) = nanstd(trials_block2_PETH(jj,:)/sqrt(numel(trials_block2_PETH(jj,:))-sum(isnan(trials_block2_PETH(jj,:)))));
        end
    end
    
    switch sesh.block_order
        case 1
            PETH.Trial.MEAN.rew_trials_light = avg_rew_trials_block1;
            PETH.Trial.MEAN.unrew_trials_light = avg_unrew_trials_block1;
            PETH.Trial.MEAN.rew_trials_sound = avg_rew_trials_block2;
            PETH.Trial.MEAN.unrew_trials_sound = avg_unrew_trials_block2;
            PETH.Trial.MEAN.trials_light_PETH = avg_trials_block1_PETH;
            PETH.Trial.MEAN.trials_sound_PETH = avg_trials_block2_PETH;
                
            PETH.Trial.SEM.rew_trials_light = SEM_rew_trials_block1;
            PETH.Trial.SEM.unrew_trials_light = SEM_unrew_trials_block1;
            PETH.Trial.SEM.rew_trials_sound = SEM_rew_trials_block2;
            PETH.Trial.SEM.unrew_trials_sound = SEM_unrew_trials_block2;
            PETH.Trial.SEM.trials_light_PETH = SEM_trials_block1_PETH; 
            PETH.Trial.SEM.trials_sound_PETH = SEM_trials_block2_PETH;
            
            PETH.Trial.ALL.rew_trials_light = rew_trials_block1;
            PETH.Trial.ALL.unrew_trials_light = unrew_trials_block1;
            PETH.Trial.ALL.rew_trials_sound = rew_trials_block2;
            PETH.Trial.ALL.unrew_trials_sound = unrew_trials_block2;
            PETH.Trial.ALL.trials_light_PETH = trials_block1_PETH;
            PETH.Trial.ALL.trials_sound_PETH = trials_block2_PETH;
            
        case 2
            PETH.Trial.MEAN.rew_trials_light = avg_rew_trials_block2;
            PETH.Trial.MEAN.unrew_trials_light = avg_unrew_trials_block2;
            PETH.Trial.MEAN.rew_trials_sound = avg_rew_trials_block1;
            PETH.Trial.MEAN.unrew_trials_sound = avg_unrew_trials_block1;
            PETH.Trial.MEAN.trials_light_PETH = avg_trials_block2_PETH;
            PETH.Trial.MEAN.trials_sound_PETH = avg_trials_block1_PETH;
            
            PETH.Trial.SEM.rew_trials_light = SEM_rew_trials_block2;
            PETH.Trial.SEM.unrew_trials_light = SEM_unrew_trials_block2;
            PETH.Trial.SEM.rew_trials_sound = SEM_rew_trials_block1;
            PETH.Trial.SEM.unrew_trials_sound = SEM_unrew_trials_block1;
            PETH.Trial.SEM.trials_light_PETH = SEM_trials_block2_PETH; 
            PETH.Trial.SEM.trials_sound_PETH = SEM_trials_block1_PETH;
            
            PETH.Trial.ALL.rew_trials_light = rew_trials_block2;
            PETH.Trial.ALL.unrew_trials_light = unrew_trials_block2;
            PETH.Trial.ALL.rew_trials_sound = rew_trials_block1;
            PETH.Trial.ALL.unrew_trials_sound = unrew_trials_block1;
            PETH.Trial.ALL.trials_light_PETH = trials_block2_PETH;
            PETH.Trial.ALL.trials_sound_PETH = trials_block1_PETH;
            
    end

    PETH.Trial.MEAN.trials_PETH = avg_trials_PETH;
    PETH.Trial.MEAN.trials_rew_PETH = avg_trials_rew_PETH;
    PETH.Trial.MEAN.trials_unrew_PETH = avg_trials_unrew_PETH;
    
    PETH.Trial.SEM.trials_PETH = SEM_trials_PETH;
    PETH.Trial.SEM.trials_rew_PETH = SEM_trials_rew_PETH;
    PETH.Trial.SEM.trials_unrew_PETH = SEM_trials_unrew_PETH;
    
    PETH.Trial.ALL.trials_PETH = trials_PETH;
    PETH.Trial.ALL.trials_rew_PETH = trials_rew_PETH;
    PETH.Trial.ALL.trials_unrew_PETH = trials_unrew_PETH;
    
end

%% trials separated by arm location
if sesh.PETH.Arm == 1
    
    temp = -5000;
    photosensor1_count = 1;
    photosensor2_count = 1;
    photosensor3_count = 1;
    photosensor4_count = 1;
    photosensor1_b1_count = 1;
    photosensor2_b1_count = 1;
    photosensor3_b1_count = 1;
    photosensor4_b1_count = 1;
    photosensor1_b2_count = 1;
    photosensor2_b2_count = 1;
    photosensor3_b2_count = 1;
    photosensor4_b2_count = 1;
    for ik = 1:length(meta.TrialInfo_block1.trialT)
        switch meta.TrialInfo_block1.photosensorID(ik) %add to count of location of trial
            case 1
                for jk = 1:15001
                    if ik == 1
                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        photosensor1_block1(jk,photosensor1_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    photosensor1(jk,photosensor1_count) = NaN;
                                    photosensor1_block1(jk,photosensor1_b1_count) = NaN;
                                else
                                    photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    photosensor1_block1(jk,photosensor1_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                photosensor1(jk,photosensor1_count) = NaN;
                                photosensor1_block1(jk,photosensor1_b1_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            photosensor1(jk:15001,photosensor1_count) = NaN;
                            photosensor1_block1(jk:15001,photosensor1_b1_count) = NaN;
                            break
                        else
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                             photosensor1_block1(jk,photosensor1_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor1_count = photosensor1_count + 1;
                photosensor1_b1_count = photosensor1_b1_count + 1;
            case 2
                for jk = 1:15001
                    if ik == 1
                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        photosensor2_block1(jk,photosensor2_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    photosensor2(jk,photosensor2_count) = NaN;
                                    photosensor2_block1(jk,photosensor2_b1_count) = NaN;
                                else
                                    photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    photosensor2_block1(jk,photosensor2_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                photosensor2(jk,photosensor2_count) = NaN;
                                photosensor2_block1(jk,photosensor2_b1_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            photosensor2(jk:15001,photosensor2_count) = NaN;
                            photosensor2_block1(jk:15001,photosensor2_b1_count) = NaN;
                            break
                        else
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            photosensor2_block1(jk,photosensor2_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor2_count = photosensor2_count + 1;
                photosensor2_b1_count = photosensor2_b1_count + 1;
            case 3
                for jk = 1:15001
                    if ik == 1
                        photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        photosensor3_block1(jk,photosensor3_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    photosensor3(jk,photosensor3_count) = NaN;
                                    photosensor3_block1(jk,photosensor3_b1_count) = NaN;
                                else
                                    photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    photosensor3_block1(jk,photosensor3_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                photosensor3(jk,photosensor3_count) = NaN;
                                photosensor3_block1(jk,photosensor3_b1_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            photosensor3(jk:15001,photosensor3_count) = NaN;
                            photosensor3_block1(jk:15001,photosensor3_b1_count) = NaN;
                            break
                        else
                            photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            photosensor3_block1(jk,photosensor3_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor3_count = photosensor3_count + 1;
                photosensor3_b1_count = photosensor3_b1_count + 1;
            case 4
                for jk = 1:15001
                    if ik == 1
                        photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        photosensor4_block1(jk,photosensor4_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    photosensor4(jk,photosensor4_count) = NaN;
                                    photosensor4_block1(jk,photosensor4_b1_count) = NaN;
                                else
                                    photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    photosensor4_block1(jk,photosensor4_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                photosensor4(jk,photosensor4_count) = NaN;
                                photosensor4_block1(jk,photosensor4_b1_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            photosensor4(jk:15001,photosensor4_count) = NaN;
                            photosensor4_block1(jk:15001,photosensor4_b1_count) = NaN;
                            break
                        else
                            photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                            photosensor4_block1(jk,photosensor4_b1_count) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor4_count = photosensor4_count + 1;
                photosensor4_b1_count = photosensor4_b1_count + 1;
        end
        temp = -5000;
    end
    
    for ip = 1:length(meta.TrialInfo_block2.trialT)
        switch meta.TrialInfo_block2.photosensorID(ip) %add to count of location of trial
            case 1
                for jk = 1:15001
                    if ip == 1
                        photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        photosensor1_block2(jk,photosensor1_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    photosensor1(jk,photosensor1_count) = NaN;
                                    photosensor1_block2(jk,photosensor1_b2_count) = NaN;
                                else
                                    photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    photosensor1_block2(jk,photosensor1_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                photosensor1(jk,photosensor1_count) = NaN;
                                photosensor1_block2(jk,photosensor1_b2_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            photosensor1(jk:15001,photosensor1_count) = NaN;
                            photosensor1_block2(jk:15001,photosensor1_b2_count) = NaN;
                            break
                        else
                            photosensor1(jk,photosensor1_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            photosensor1_block2(jk,photosensor1_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor1_count = photosensor1_count + 1;
                photosensor1_b2_count = photosensor1_b2_count + 1;
            case 2
                for jk = 1:15001
                    if ip == 1
                        photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        photosensor2_block2(jk,photosensor2_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    photosensor2(jk,photosensor2_count) = NaN;
                                    photosensor2_block2(jk,photosensor2_b2_count) = NaN;
                                else
                                    photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                     photosensor2_block2(jk,photosensor2_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                photosensor2(jk,photosensor2_count) = NaN;
                                photosensor2_block2(jk,photosensor2_b2_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            photosensor2(jk:15001,photosensor2_count) = NaN;
                            photosensor2_block2(jk:15001,photosensor2_b2_count) = NaN;
                            break
                        else
                            photosensor2(jk,photosensor2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                             photosensor2_block2(jk,photosensor2_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor2_count = photosensor2_count + 1;
                photosensor2_b2_count = photosensor2_b2_count + 1;
            case 3
                for jk = 1:15001
                    if ip == 1
                        photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        photosensor3_block2(jk,photosensor3_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    photosensor3(jk,photosensor3_count) = NaN;
                                    photosensor3_block2(jk,photosensor3_b2_count) = NaN;
                                else
                                    photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    photosensor3_block2(jk,photosensor3_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                photosensor3(jk,photosensor3_count) = NaN;
                                photosensor3_block2(jk,photosensor3_b2_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            photosensor3(jk:15001,photosensor3_count) = NaN;
                            photosensor3_block2(jk:15001,photosensor3_b2_count) = NaN;
                            break
                        else
                            photosensor3(jk,photosensor3_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            photosensor3_block2(jk,photosensor3_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor3_count = photosensor3_count + 1;
                photosensor3_b2_count = photosensor3_b2_count + 1;
            case 4
                for jk = 1:15001
                    if ip == 1
                        photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        photosensor4_block2(jk,photosensor4_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    photosensor4(jk,photosensor4_count) = NaN;
                                    photosensor4_block2(jk,photosensor4_b2_count) = NaN;
                                else
                                    photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    photosensor4_block2(jk,photosensor4_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                photosensor4(jk,photosensor4_count) = NaN;
                                photosensor4_block2(jk,photosensor4_b2_count) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            photosensor4(jk:15001,photosensor4_count) = NaN;
                            photosensor4_block2(jk:15001,photosensor4_b2_count) = NaN;
                            break
                        else
                            photosensor4(jk,photosensor4_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            photosensor4_block2(jk,photosensor4_b2_count) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                photosensor4_count = photosensor4_count + 1;
                photosensor4_b2_count = photosensor4_b2_count + 1;
        end
        temp = -5000;
    end
    
    for jj = 1:15001 %generate averaged responses
        if jj <= length(photosensor1)
            PETH.Arm.MEAN.photosensor1(jj,1) = nanmean(photosensor1(jj,:));
            PETH.Arm.SEM.photosensor1(jj,1) = nanstd(photosensor1(jj,:)/sqrt(numel(photosensor1(jj,:))-sum(isnan(photosensor1(jj,:)))));
        end
        if jj <= length(photosensor2)
            PETH.Arm.MEAN.photosensor2(jj,1) = nanmean(photosensor2(jj,:));
            PETH.Arm.SEM.photosensor2(jj,1) = nanstd(photosensor2(jj,:)/sqrt(numel(photosensor2(jj,:))-sum(isnan(photosensor2(jj,:)))));
        end
        if jj <= length(photosensor3)
            PETH.Arm.MEAN.photosensor3(jj,1) = nanmean(photosensor3(jj,:));
            PETH.Arm.SEM.photosensor3(jj,1) = nanstd(photosensor3(jj,:)/sqrt(numel(photosensor3(jj,:))-sum(isnan(photosensor3(jj,:)))));
        end
        if jj <= length(photosensor4)
            PETH.Arm.MEAN.photosensor4(jj,1) = nanmean(photosensor4(jj,:));
            PETH.Arm.SEM.photosensor4(jj,1) = nanstd(photosensor4(jj,:)/sqrt(numel(photosensor4(jj,:))-sum(isnan(photosensor4(jj,:)))));
        end
       if jj <= length(photosensor1_block1)
            avg_photosensor1_block1(jj,1) = nanmean(photosensor1_block1(jj,:));
            SEM_photosensor1_block1(jj,1) = nanstd(photosensor1_block1(jj,:)/sqrt(numel(photosensor1_block1(jj,:))-sum(isnan(photosensor1_block1(jj,:)))));
        end
        if jj <= length(photosensor2_block1)
            avg_photosensor2_block1(jj,1) = nanmean(photosensor2_block1(jj,:));
            SEM_photosensor2_block1(jj,1) = nanstd(photosensor2_block1(jj,:)/sqrt(numel(photosensor2_block1(jj,:))-sum(isnan(photosensor2_block1(jj,:)))));
        end
        if jj <= length(photosensor3_block1)
            avg_photosensor3_block1(jj,1) = nanmean(photosensor3_block1(jj,:));
            SEM_photosensor3_block1(jj,1) = nanstd(photosensor3_block1(jj,:)/sqrt(numel(photosensor3_block1(jj,:))-sum(isnan(photosensor3_block1(jj,:)))));
        end
        if jj <= length(photosensor4_block1)
            avg_photosensor4_block1(jj,1) = nanmean(photosensor4_block1(jj,:));
            SEM_photosensor4_block1(jj,1) = nanstd(photosensor4_block1(jj,:)/sqrt(numel(photosensor4_block1(jj,:))-sum(isnan(photosensor4_block1(jj,:)))));
        end
if jj <= length(photosensor1_block2)
            avg_photosensor1_block2(jj,1) = nanmean(photosensor1_block2(jj,:));
            SEM_photosensor1_block2(jj,1) = nanstd(photosensor1_block2(jj,:)/sqrt(numel(photosensor1_block2(jj,:))-sum(isnan(photosensor1_block2(jj,:)))));
        end
        if jj <= length(photosensor2_block2)
            avg_photosensor2_block2(jj,1) = nanmean(photosensor2_block2(jj,:));
            SEM_photosensor2_block2(jj,1) = nanstd(photosensor2_block2(jj,:)/sqrt(numel(photosensor2_block2(jj,:))-sum(isnan(photosensor2_block2(jj,:)))));
        end
        if jj <= length(photosensor3_block2)
            avg_photosensor3_block2(jj,1) = nanmean(photosensor3_block2(jj,:));
            SEM_photosensor3_block2(jj,1) = nanstd(photosensor3_block2(jj,:)/sqrt(numel(photosensor3_block2(jj,:))-sum(isnan(photosensor3_block2(jj,:)))));
        end
        if jj <= length(photosensor4_block2)
            avg_photosensor4_block2(jj,1) = nanmean(photosensor4_block2(jj,:));
            SEM_photosensor4_block2(jj,1) = nanstd(photosensor4_block2(jj,:)/sqrt(numel(photosensor4_block2(jj,:))-sum(isnan(photosensor4_block2(jj,:)))));
        end
    end
    
    switch sesh.block_order
        case 1
            PETH.Arm.MEAN.photosensor1_light = avg_photosensor1_block1;
            PETH.Arm.MEAN.photosensor2_light = avg_photosensor2_block1;
            PETH.Arm.MEAN.photosensor3_light = avg_photosensor3_block1;
            PETH.Arm.MEAN.photosensor4_light = avg_photosensor4_block1;
            PETH.Arm.MEAN.photosensor1_sound = avg_photosensor1_block2;
            PETH.Arm.MEAN.photosensor2_sound = avg_photosensor2_block2;
            PETH.Arm.MEAN.photosensor3_sound = avg_photosensor3_block2;
            PETH.Arm.MEAN.photosensor4_sound = avg_photosensor4_block2;
                
            PETH.Arm.SEM.photosensor1_light = SEM_photosensor1_block1;
            PETH.Arm.SEM.photosensor2_light = SEM_photosensor2_block1;
            PETH.Arm.SEM.photosensor3_light = SEM_photosensor3_block1;
            PETH.Arm.SEM.photosensor4_light = SEM_photosensor4_block1;
            PETH.Arm.SEM.photosensor1_sound = SEM_photosensor1_block2;
            PETH.Arm.SEM.photosensor2_sound = SEM_photosensor2_block2;
            PETH.Arm.SEM.photosensor3_sound = SEM_photosensor3_block2;
            PETH.Arm.SEM.photosensor4_sound = SEM_photosensor4_block2;           
            
        case 2
            PETH.Arm.MEAN.photosensor1_light = avg_photosensor1_block2;
            PETH.Arm.MEAN.photosensor2_light = avg_photosensor2_block2;
            PETH.Arm.MEAN.photosensor3_light = avg_photosensor3_block2;
            PETH.Arm.MEAN.photosensor4_light = avg_photosensor4_block2;
            PETH.Arm.MEAN.photosensor1_sound = avg_photosensor1_block1;
            PETH.Arm.MEAN.photosensor2_sound = avg_photosensor2_block1;
            PETH.Arm.MEAN.photosensor3_sound = avg_photosensor3_block1;
            PETH.Arm.MEAN.photosensor4_sound = avg_photosensor4_block1;
                
            PETH.Arm.SEM.photosensor1_light = SEM_photosensor1_block2;
            PETH.Arm.SEM.photosensor2_light = SEM_photosensor2_block2;
            PETH.Arm.SEM.photosensor3_light = SEM_photosensor3_block2;
            PETH.Arm.SEM.photosensor4_light = SEM_photosensor4_block2;
            PETH.Arm.SEM.photosensor1_sound = SEM_photosensor1_block1;
            PETH.Arm.SEM.photosensor2_sound = SEM_photosensor2_block1;
            PETH.Arm.SEM.photosensor3_sound = SEM_photosensor3_block1;
            PETH.Arm.SEM.photosensor4_sound = SEM_photosensor4_block1; 
            
    end
    
    PETH.Arm.ALL.photosensor1 = photosensor1;
    PETH.Arm.ALL.photosensor2 = photosensor2;
    PETH.Arm.ALL.photosensor3 = photosensor3;
    PETH.Arm.ALL.photosensor4 = photosensor4;
    
end

%% PETHs aligned to cue-onset for trials where the rat NOSEPOKED
if sesh.PETH.Trial_np == 1
    
    temp = -5000;
    rew_trial_num_block1 = 1;
    unrew_trial_num_block1 = 1;
    rew_trial_num_block2 = 1;
    unrew_trial_num_block2 = 1;
    trial_num = 1;
    trial_num2 = 1;
    for ik = 1:length(meta.TrialInfo_block1.trialT)
        if meta.TrialInfo_block1.nosepokeT(ik) ~= 0
        switch meta.TrialInfo_block1.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                for jk = 1:15001
                    if ik == 1
                        rew_trials4np_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    rew_trials4np_block1(jk,rew_trial_num_block1) = NaN;
                                    trials4np_PETH(jk,trial_num) = NaN;
                                    trials4np_block1_PETH(jk,trial_num) = NaN;
                                else
                                    rew_trials4np_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                rew_trials4np_block1(jk,rew_trial_num_block1) = NaN;
                                trials4np_PETH(jk,trial_num) = NaN;
                                trials4np_block1_PETH(jk,trial_num) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            rew_trials4np_block1(jk:15001,rew_trial_num_block1) = NaN;
                            trials4np_PETH(jk:15001,trial_num) = NaN;
                            trials4np_block1_PETH(jk:15001,trial_num) = NaN;
                            break
                        else
                            rew_trials4np_block1(jk,rew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                rew_trial_num_block1 = rew_trial_num_block1 + 1;
                trial_num = trial_num + 1;
            case 0
                for jk = 1:15001
                    if ik == 1
                        unrew_trials4np_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                            if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                    unrew_trials4np_block1(jk,unrew_trial_num_block1) = NaN;
                                    trials4np_PETH(jk,trial_num) = NaN;
                                    trials4np_block1_PETH(jk,trial_num) = NaN;
                                else
                                    unrew_trials4np_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                                end
                            else
                                unrew_trials4np_block1(jk,unrew_trial_num_block1) = NaN;
                                trials4np_PETH(jk,trial_num) = NaN;
                                trials4np_block1_PETH(jk,trial_num) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            unrew_trials4np_block1(jk:15001,unrew_trial_num_block1) = NaN;
                            trials4np_PETH(jk:15001,trial_num) = NaN;
                            trials4np_block1_PETH(jk:15001,trial_num) = NaN;
                            break
                        else
                            unrew_trials4np_block1(jk,unrew_trial_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            trials4np_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
                trial_num = trial_num + 1;
        end
        
        temp = -5000;
        end
    end
    
    for ip = 1:length(meta.TrialInfo_block2.trialT)
        if meta.TrialInfo_block2.nosepokeT(ip) ~= 0
        switch meta.TrialInfo_block2.rewarded(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                for jk = 1:15001
                    if ip == 1
                        rew_trials4np_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    rew_trials4np_block2(jk,rew_trial_num_block2) = NaN;
                                    trials4np_PETH(jk,trial_num) = NaN;
                                    trials4np_block2_PETH(jk,trial_num2) = NaN;
                                else
                                    rew_trials4np_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                rew_trials4np_block2(jk,rew_trial_num_block2) = NaN;
                                trials4np_PETH(jk,trial_num) = NaN;
                                trials4np_block2_PETH(jk,trial_num2) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            rew_trials4np_block2(jk:15001,rew_trial_num_block2) = NaN;
                            trials4np_PETH(jk:15001,trial_num) = NaN;
                            trials4np_block2_PETH(jk:15001,trial_num2) = NaN;
                            break
                        else
                            rew_trials4np_block2(jk,rew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                rew_trial_num_block2 = rew_trial_num_block2 + 1;
                trial_num = trial_num + 1;
                trial_num2 = trial_num2 + 1;
            case 0
                for jk = 1:15001
                    if ip == 1
                        unrew_trials4np_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                    else
                        if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                            if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                    unrew_trials4np_block2(jk,unrew_trial_num_block2) = NaN;
                                    trials4np_PETH(jk,trial_num) = NaN;
                                    trials4np_block2_PETH(jk,trial_num2) = NaN;
                                else
                                    unrew_trials4np_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                end
                            else
                                unrew_trials4np_block2(jk,unrew_trial_num_block2) = NaN;
                                trials4np_PETH(jk,trial_num) = NaN;
                                trials4np_block2_PETH(jk,trial_num2) = NaN;
                            end
                        elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                            unrew_trials4np_block2(jk:15001,unrew_trial_num_block2) = NaN;
                            trials4np_PETH(jk:15001,trial_num) = NaN;
                            trials4np_block2_PETH(jk:15001,trial_num2) = NaN;
                            break
                        else
                            unrew_trials4np_block2(jk,unrew_trial_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials4np_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            trials4np_block2_PETH(jk,trial_num2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        end
                    end
                    temp = temp + 1;
                end
                unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
                trial_num = trial_num + 1;
                trial_num2 = trial_num2 + 1;
        end
        temp = -5000;
        end
    end
    
    trials4np_rew_PETH = cat(2,rew_trials4np_block1,rew_trials4np_block2);
    trials4np_unrew_PETH = cat(2,unrew_trials4np_block1,unrew_trials4np_block2);
    
    for jj = 1:15001 %generate averaged responses
        if jj <= length(rew_trials4np_block1)
            avg_rew_trials4np_block1(jj,1) = nanmean(rew_trials4np_block1(jj,:));
            SEM_rew_trials4np_block1(jj,1) = nanstd(rew_trials4np_block1(jj,:)/sqrt(numel(rew_trials4np_block1(jj,:))-sum(isnan(rew_trials4np_block1(jj,:)))));
        end
        if jj <= length(unrew_trials4np_block1)
            avg_unrew_trials4np_block1(jj,1) = nanmean(unrew_trials4np_block1(jj,:));
            SEM_unrew_trials4np_block1(jj,1) = nanstd(unrew_trials4np_block1(jj,:)/sqrt(numel(unrew_trials4np_block1(jj,:))-sum(isnan(unrew_trials4np_block1(jj,:)))));
        end
        if jj <= length(rew_trials4np_block2)
            avg_rew_trials4np_block2(jj,1) = nanmean(rew_trials4np_block2(jj,:));
            SEM_rew_trials4np_block2(jj,1) = nanstd(rew_trials4np_block2(jj,:)/sqrt(numel(rew_trials4np_block2(jj,:))-sum(isnan(rew_trials4np_block2(jj,:)))));
        end
        if jj <= length(unrew_trials4np_block2)
            avg_unrew_trials4np_block2(jj,1) = nanmean(unrew_trials4np_block2(jj,:));
            SEM_unrew_trials4np_block2(jj,1) = nanstd(unrew_trials4np_block2(jj,:)/sqrt(numel(unrew_trials4np_block2(jj,:))-sum(isnan(unrew_trials4np_block2(jj,:)))));
        end
        
        if jj <= length(trials4np_rew_PETH)
            avg_trials4np_rew_PETH(jj,1) = nanmean(trials4np_rew_PETH(jj,:));
            SEM_trials4np_rew_PETH(jj,1) = nanstd(trials4np_rew_PETH(jj,:)/sqrt(numel(trials4np_rew_PETH(jj,:))-sum(isnan(trials4np_rew_PETH(jj,:)))));
        end
        
        if jj <= length(trials4np_unrew_PETH)
            avg_trials4np_unrew_PETH(jj,1) = nanmean(trials4np_unrew_PETH(jj,:));
            SEM_trials4np_unrew_PETH(jj,1) = nanstd(trials4np_unrew_PETH(jj,:)/sqrt(numel(trials4np_unrew_PETH(jj,:))-sum(isnan(trials4np_unrew_PETH(jj,:)))));
        end
        
        if jj <= length(trials4np_PETH)
            avg_trials4np_PETH(jj,1) = nanmean(trials4np_PETH(jj,:));
            SEM_trials4np_PETH(jj,1) = nanstd(trials4np_PETH(jj,:)/sqrt(numel(trials4np_PETH(jj,:))-sum(isnan(trials4np_PETH(jj,:)))));
        end
        
        if jj <= length(trials4np_block1_PETH)
            avg_trials4np_block1_PETH(jj,1) = nanmean(trials4np_block1_PETH(jj,:));
            SEM_trials4np_block1_PETH(jj,1) = nanstd(trials4np_block1_PETH(jj,:)/sqrt(numel(trials4np_block1_PETH(jj,:))-sum(isnan(trials4np_block1_PETH(jj,:)))));
        end
        
        if jj <= length(trials4np_block2_PETH)
            avg_trials4np_block2_PETH(jj,1) = nanmean(trials4np_block2_PETH(jj,:));
            SEM_trials4np_block2_PETH(jj,1) = nanstd(trials4np_block2_PETH(jj,:)/sqrt(numel(trials4np_block2_PETH(jj,:))-sum(isnan(trials4np_block2_PETH(jj,:)))));
        end
    end
    
    switch sesh.block_order
        case 1
            PETH.Trial4np.MEAN.rew_trials4np_light = avg_rew_trials4np_block1;
            PETH.Trial4np.MEAN.unrew_trials4np_light = avg_unrew_trials4np_block1;
            PETH.Trial4np.MEAN.rew_trials4np_sound = avg_rew_trials4np_block2;
            PETH.Trial4np.MEAN.unrew_trials4np_sound = avg_unrew_trials4np_block2;
            PETH.Trial4np.MEAN.trials4np_light_PETH = avg_trials4np_block1_PETH;
            PETH.Trial4np.MEAN.trials4np_sound_PETH = avg_trials4np_block2_PETH;
                
            PETH.Trial4np.SEM.rew_trials4np_light = SEM_rew_trials4np_block1;
            PETH.Trial4np.SEM.unrew_trials4np_light = SEM_unrew_trials4np_block1;
            PETH.Trial4np.SEM.rew_trials4np_sound = SEM_rew_trials4np_block2;
            PETH.Trial4np.SEM.unrew_trials4np_sound = SEM_unrew_trials4np_block2;
            PETH.Trial4np.SEM.trials4np_light_PETH = SEM_trials4np_block1_PETH; 
            PETH.Trial4np.SEM.trials4np_sound_PETH = SEM_trials4np_block2_PETH;
            
            PETH.Trial4np.ALL.rew_trials4np_light = rew_trials4np_block1;
            PETH.Trial4np.ALL.unrew_trials4np_light = unrew_trials4np_block1;
            PETH.Trial4np.ALL.rew_trials4np_sound = rew_trials4np_block2;
            PETH.Trial4np.ALL.unrew_trials4np_sound = unrew_trials4np_block2;
            PETH.Trial4np.ALL.trials4np_light_PETH = trials4np_block1_PETH;
            PETH.Trial4np.ALL.trials4np_sound_PETH = trials4np_block2_PETH;
            
        case 2
            PETH.Trial4np.MEAN.rew_trials4np_light = avg_rew_trials4np_block2;
            PETH.Trial4np.MEAN.unrew_trials4np_light = avg_unrew_trials4np_block2;
            PETH.Trial4np.MEAN.rew_trials4np_sound = avg_rew_trials4np_block1;
            PETH.Trial4np.MEAN.unrew_trials4np_sound = avg_unrew_trials4np_block1;
            PETH.Trial4np.MEAN.trials4np_light_PETH = avg_trials4np_block2_PETH;
            PETH.Trial4np.MEAN.trials4np_sound_PETH = avg_trials4np_block1_PETH;
            
            PETH.Trial4np.SEM.rew_trials4np_light = SEM_rew_trials4np_block2;
            PETH.Trial4np.SEM.unrew_trials4np_light = SEM_unrew_trials4np_block2;
            PETH.Trial4np.SEM.rew_trials4np_sound = SEM_rew_trials4np_block1;
            PETH.Trial4np.SEM.unrew_trials4np_sound = SEM_unrew_trials4np_block1;
            PETH.Trial4np.SEM.trials4np_light_PETH = SEM_trials4np_block2_PETH; 
            PETH.Trial4np.SEM.trials4np_sound_PETH = SEM_trials4np_block1_PETH;
            
            PETH.Trial4np.ALL.rew_trials4np_light = rew_trials4np_block2;
            PETH.Trial4np.ALL.unrew_trials4np_light = unrew_trials4np_block2;
            PETH.Trial4np.ALL.rew_trials4np_sound = rew_trials4np_block1;
            PETH.Trial4np.ALL.unrew_trials4np_sound = unrew_trials4np_block1;
            PETH.Trial4np.ALL.trials4np_light_PETH = trials4np_block2_PETH;
            PETH.Trial4np.ALL.trials4np_sound_PETH = trials4np_block1_PETH;
            
    end

    PETH.Trial4np.MEAN.trials4np_PETH = avg_trials4np_PETH;
    PETH.Trial4np.MEAN.trials4np_rew_PETH = avg_trials4np_rew_PETH;
    PETH.Trial4np.MEAN.trials4np_unrew_PETH = avg_trials4np_unrew_PETH;
    
    PETH.Trial4np.SEM.trials4np_PETH = SEM_trials4np_PETH;
    PETH.Trial4np.SEM.trials4np_rew_PETH = SEM_trials4np_rew_PETH;
    PETH.Trial4np.SEM.trials4np_unrew_PETH = SEM_trials4np_unrew_PETH;
    
    PETH.Trial4np.ALL.trials4np_PETH = trials4np_PETH;
    PETH.Trial4np.ALL.trials4np_rew_PETH = trials4np_rew_PETH;
    PETH.Trial4np.ALL.trials4np_unrew_PETH = trials4np_unrew_PETH;
    
end

%% unrewarded approach vs skip trials
if sesh.PETH.Approach == 1
    
    temp = -5000;
    skip_num_block1 = 1;
    app_num_block1 = 1;
    skip_num_block2 = 1;
    app_num_block2 = 1;
    for ik = 1:length(meta.TrialInfo_block1.nosepokeT)
        if meta.TrialInfo_block1.rewarded(ik) == 0
            switch meta.TrialInfo_block1.nosepokeID(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
                case 0
                    for jk = 1:15001
                        if ik == 1
                            unrew_skip_block1(jk,skip_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                                if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                        unrew_skip_block1(jk,skip_num_block1) = NaN;
                                    else
                                        unrew_skip_block1(jk,skip_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    unrew_skip_block1(jk,skip_num_block1) = NaN;
                                end
                            elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                                unrew_skip_block1(jk:15001,skip_num_block1) = NaN;
                                break
                            else
                                unrew_skip_block1(jk,skip_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    skip_num_block1 = skip_num_block1 + 1;
                case {5,6,7,8}
                    for jk = 1:15001
                        if ik == 1
                            unrew_app_block1(jk,app_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                        else
                            if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                                if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                                    if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                        unrew_app_block1(jk,app_num_block1) = NaN;
                                    else
                                        unrew_app_block1(jk,app_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);
                                    end
                                else
                                    unrew_app_block1(jk,app_num_block1) = NaN;
                                end                       
                            elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                            unrew_app_block1(jk:15001,app_num_block1) = NaN;
                            break
                            else
                           unrew_app_block1(jk,app_num_block1) = gau_sdf(dataPoint.Trials(ik) + temp);     
                        end
                        end
                        temp = temp + 1;
                    end
                    app_num_block1 = app_num_block1 + 1;
            end
            temp = -5000;
        end
    end
    
    for ip = 1:length(meta.TrialInfo_block2.nosepokeT)
        if meta.TrialInfo_block2.rewarded(ip) == 0
            switch meta.TrialInfo_block2.nosepokeID(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
                case 0
                    for jk = 1:15001
                        if ip == 1
                            unrew_skip_block2(jk,skip_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        else
                            if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                                if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                        unrew_skip_block2(jk,skip_num_block2) = NaN;
                                    else
                                        unrew_skip_block2(jk,skip_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    end
                                else
                                    unrew_skip_block2(jk,skip_num_block2) = NaN;
                                end
                            elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                                unrew_skip_block2(jk:15001,skip_num_block2) = NaN;
                                break
                            else
                                unrew_skip_block2(jk,skip_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    skip_num_block2 = skip_num_block2 + 1;
                case {5,6,7,8}
                    for jk = 1:15001
                        if ip == 1
                            unrew_app_block2(jk,app_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                        else
                            if temp / -1000 > meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1)
                                if meta.TrialInfo_block2.unnosepoke_to_trialT(ip-1) == 0
                                    if temp / -1000 > meta.TrialInfo_block2.summary(ip-1,9)
                                        unrew_app_block2(jk,app_num_block2) = NaN;
                                    else
                                        unrew_app_block2(jk,app_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                                    end
                                else
                                    unrew_app_block2(jk,app_num_block2) = NaN;
                                end
                            elseif temp / 1000 > meta.TrialInfo_block2.summary(ip,9)
                                unrew_app_block2(jk:15001,app_num_block2) = NaN;
                                break
                            else
                                unrew_app_block2(jk,app_num_block2) = gau_sdf(dataPoint.Trials(length(meta.TrialInfo_block1.trialT) + ip) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    app_num_block2 = app_num_block2 + 1;
            end
            temp = -5000;
        end
    end
    
    for jj = 1:15001 %generate averaged responses
        if jj <= length(unrew_skip_block1)
            avg_unrew_skip_block1(jj,1) = nanmean(unrew_skip_block1(jj,:));
            SEM_unrew_skip_block1(jj,1) = nanstd(unrew_skip_block1(jj,:)/sqrt(numel(unrew_skip_block1(jj,:))-sum(isnan(unrew_skip_block1(jj,:)))));
        end
        if jj <= length(unrew_app_block1)
            avg_unrew_app_block1(jj,1) = nanmean(unrew_app_block1(jj,:));
            SEM_unrew_app_block1(jj,1) = nanstd(unrew_app_block1(jj,:)/sqrt(numel(unrew_app_block1(jj,:))-sum(isnan(unrew_app_block1(jj,:)))));
        end
        if jj <= length(unrew_skip_block2)
            avg_unrew_skip_block2(jj,1) = nanmean(unrew_skip_block2(jj,:));
            SEM_unrew_skip_block2(jj,1) = nanstd(unrew_skip_block2(jj,:)/sqrt(numel(unrew_skip_block2(jj,:))-sum(isnan(unrew_skip_block2(jj,:)))));
        end
        if jj <= length(unrew_app_block2)
            avg_unrew_app_block2(jj,1) = nanmean(unrew_app_block2(jj,:));
            SEM_unrew_app_block2(jj,1) = nanstd(unrew_app_block2(jj,:)/sqrt(numel(unrew_app_block2(jj,:))-sum(isnan(unrew_app_block2(jj,:)))));
        end
    end
    
    switch sesh.block_order
        case 1
            PETH.Approach.MEAN.unrew_skip_light = avg_unrew_skip_block1;
            PETH.Approach.MEAN.unrew_app_light = avg_unrew_app_block1;
            PETH.Approach.MEAN.unrew_skip_sound = avg_unrew_skip_block2;
            PETH.Approach.MEAN.unrew_app_sound = avg_unrew_app_block2;
            
            PETH.Approach.SEM.unrew_skip_light = SEM_unrew_skip_block1;
            PETH.Approach.SEM.unrew_app_light = SEM_unrew_app_block1;
            PETH.Approach.SEM.unrew_skip_sound = SEM_unrew_skip_block2;
            PETH.Approach.SEM.unrew_app_sound = SEM_unrew_app_block2;
        case 2
            PETH.Approach.MEAN.unrew_skip_light = avg_unrew_skip_block2;
            PETH.Approach.MEAN.unrew_app_light = avg_unrew_app_block2;
            PETH.Approach.MEAN.unrew_skip_sound = avg_unrew_skip_block1;
            PETH.Approach.MEAN.unrew_app_sound = avg_unrew_app_block1;
            
            PETH.Approach.SEM.unrew_skip_light = SEM_unrew_skip_block2;
            PETH.Approach.SEM.unrew_app_light = SEM_unrew_app_block2;
            PETH.Approach.SEM.unrew_skip_sound = SEM_unrew_skip_block1;
            PETH.Approach.SEM.unrew_app_sound = SEM_unrew_app_block1;
    end
end
toc