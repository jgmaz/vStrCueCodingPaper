%%
t = [0 4000];
binsize = 0.001; % in seconds
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;

spk_count = histc(spk_t,tbin_edges);
spk_count = spk_count(1:end-1);

gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.1./binsize; % 0.1 seconds (100ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window

%%
temp = -5000;
trial_num = 1;
trial_num2 = 1;
for ik = 1:length(meta.TrialInfo_block1.trialT)
    switch meta.TrialInfo_block1.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
            for jk = 1:15001
                if ik == 1
                    trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                else
                    if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                        if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                            if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                trials_block1_PETH(jk,trial_num) = NaN;
                            else
                                trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        else
                            trials_block1_PETH(jk,trial_num) = NaN;
                        end
                    elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                        trials_block1_PETH(jk:15001,trial_num) = NaN;
                        break
                    else
                        trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    end
                end
                temp = temp + 1;
            end
            trial_num = trial_num + 1;
        case 0
            for jk = 1:15001
                if ik == 1
                    trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                else
                    if temp / -1000 > meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1)
                        if meta.TrialInfo_block1.unnosepoke_to_trialT(ik-1) == 0
                            if temp / -1000 > meta.TrialInfo_block1.summary(ik-1,9)
                                trials_block1_PETH(jk,trial_num) = NaN;
                            else
                                trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                            end
                        else
                            trials_block1_PETH(jk,trial_num) = NaN;
                        end
                    elseif temp / 1000 > meta.TrialInfo_block1.summary(ik,9)
                        trials_block1_PETH(jk:15001,trial_num) = NaN;
                        break
                    else
                        trials_block1_PETH(jk,trial_num) = gau_sdf(dataPoint.Trials(ik) + temp);
                    end
                end
                temp = temp + 1;
            end
            trial_num = trial_num + 1;
    end
    
    temp = -5000;
end

%%
for jj = 1:15001 %generate averaged responses
    if jj <= length(trials_block1_PETH)
        avg_trials_block1_PETH(jj,1) = nanmean(trials_block1_PETH(jj,:));
    end
end

PETH.Trial.MEAN.trials_light_PETH = avg_trials_block1_PETH;

%%
cfg.start = 4001;
cfg.end = 6000;

for kk = 1:length(dir('*.mat'))
    load(mat_files(kk).name);  
    [maximum_light(kk),max_location_light(kk)] = max(PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end));
end

[sortPETH.light.sorted,sortPETH.light.index_sort] = sort(max_location_light);

for jj = 1:length(dir('*.mat'))
    load(mat_files(sortPETH.light.index_sort(jj)).name);
    sortedPETH.light.unaltered(jj,:) = PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end);
end

for jj = 1:length(dir('*.mat'))
    mean_peth.light(jj,1) = mean(sortedPETH.light.unaltered(jj,:));
    std_peth.light(jj,1) = std(sortedPETH.light.unaltered(jj,:));
    sortedPETH.light.zscore(jj,:) = (sortedPETH.light.unaltered(jj,:) - mean_peth.light(jj)) / std_peth.light(jj);
end

imagesc(-1:.001:1,1:length(dir('*.mat')),sortedPETH.light.zscore);
