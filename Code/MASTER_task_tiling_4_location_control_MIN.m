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
    
    %     block_drift.block1_length(kk) = length(FRATE.Cue.Trial_firing_rate_block1);
    %     block_drift.block1_half(kk) = round(block_drift.block1_length(kk) / 2);
    %     block_drift.b1_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)));
    %     block_drift.b1_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    %     block_drift.MWU_b1(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    %
    %     block_drift.block2_length(kk) = length(FRATE.Cue.Trial_firing_rate_block2);
    %     block_drift.block2_half(kk) = round(block_drift.block2_length(kk) / 2);
    %     block_drift.b2_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)));
    %     block_drift.b2_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    %     block_drift.MWU_b2(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    %
    %     switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
    %                 case 0
    
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
                    
                end
                temp = -5000;
            end
    end
    PETH.Arm.ALL.photosensor1 = photosensor1;
    
    
    %%
    rng('shuffle')
    num_trials = round(length(PETH.Arm.ALL.photosensor1(1,:)));
    temp_idx = randperm(num_trials);
    half_1 = temp_idx(1:num_trials/2);
    half_2 = temp_idx(num_trials/2+1:end);
    
    for iBin = 1:15001 %generate averaged responses
        if iBin <= length(PETH.Arm.ALL.photosensor1)
            PETHs.trials_arm1_1st_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor1(iBin,half_1));
            PETHs.trials_arm1_2nd_half{count}(iBin,:) = nanmean(PETH.Arm.ALL.photosensor1(iBin,half_2));
        end
    end
    
    [minimum_arm1_1st(count),min_location_arm1_1st(count)] = min(PETHs.trials_arm1_1st_half{count}(cfg.start:cfg.end));
    [minimum_arm1_2nd(count),min_location_arm1_2nd(count)] = min(PETHs.trials_arm1_2nd_half{count}(cfg.start:cfg.end));
    
    count = count + 1;
    
    %         case 1
    %     end
end

[MINsortPETH.arm1_1st.sorted,MINsortPETH.arm1_1st.index_sort] = sort(min_location_arm1_1st);
[MINsortPETH.arm1_2nd.sorted,MINsortPETH.arm1_2nd.index_sort] = sort(min_location_arm1_2nd);

for jj = 1:length(mat_files_drift)
    MINsortedPETH.arm1_1st.unaltered(jj,:) = PETHs.trials_arm1_1st_half{MINsortPETH.arm1_1st.index_sort(jj)}(cfg.start:cfg.end);
    MINsortedPETH.arm1_2nd.arm1_1st(jj,:) = PETHs.trials_arm1_2nd_half{MINsortPETH.arm1_1st.index_sort(jj)}(cfg.start:cfg.end);
    MINsortedPETH.arm1_2nd.unaltered(jj,:) = PETHs.trials_arm1_2nd_half{MINsortPETH.arm1_2nd.index_sort(jj)}(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(mat_files_drift)
    mean_peth.arm1_1st(jj,1) = mean(MINsortedPETH.arm1_1st.unaltered(jj,:));
    std_peth.arm1_1st(jj,1) = std(MINsortedPETH.arm1_1st.unaltered(jj,:));
    MINsortedPETH.arm1_1st.zscore(jj,:) = (MINsortedPETH.arm1_1st.unaltered(jj,:) - mean_peth.arm1_1st(jj)) / std_peth.arm1_1st(jj);
    mean_peth.arm1_2nd(jj,1) = mean(MINsortedPETH.arm1_2nd.unaltered(jj,:));
    std_peth.arm1_2nd(jj,1) = std(MINsortedPETH.arm1_2nd.unaltered(jj,:));
    MINsortedPETH.arm1_2nd.zscore(jj,:) = (MINsortedPETH.arm1_2nd.unaltered(jj,:) - mean_peth.arm1_2nd(jj)) / std_peth.arm1_2nd(jj);
    mean_peth.arm1_2nd_2_1st(jj,1) = mean(MINsortedPETH.arm1_2nd.arm1_1st(jj,:));
    std_peth.arm1_2nd_2_1st(jj,1) = std(MINsortedPETH.arm1_2nd.arm1_1st(jj,:));
    MINsortedPETH.arm1_2nd_2_1st.zscore(jj,:) = (MINsortedPETH.arm1_2nd.arm1_1st(jj,:) - mean_peth.arm1_2nd_2_1st(jj)) / std_peth.arm1_2nd_2_1st(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','MINsortPETH','MINsortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:1,1:length(dir('*.mat')),MINsortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.arm1_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.arm1_2nd_2_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.arm1_2nd.zscore);
colorbar; caxis([-3 4]);