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
    receptacle1_count = 1;
    nosepoke_number = 1;
    
    new_v_old = strcmp(mat_overview.fname{kk}(1:4),'R060');
    switch new_v_old
        case 0
            %% old rats (R053,R056,R057)
    for ik = 1:length(metadata.TrialInfo_block1.nosepokeT)
        if metadata.TrialInfo_block1.nosepokeT(ik) ~= 0
            if metadata.TrialInfo_block1.nosepokeID(ik) == 5 %add to nosepoke counts
                    for jk = 1:15001
                        if ik == 1
                            receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block1.trial_to_nosepokeT(ik)
                                receptacle1(jk,receptacle1_count) = NaN;
                            elseif temp / 1000 > metadata.TrialInfo_block1.nosepoke_length(ik)
                                receptacle1(jk:15001,receptacle1_count) = NaN;
                                break
                            else
                                receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    receptacle1_count = receptacle1_count + 1;

            end
            temp = -5000;
            nosepoke_number = nosepoke_number + 1;
        end
    end
    
    for ip = 1:length(metadata.TrialInfo_block2.nosepokeT)
        if metadata.TrialInfo_block2.nosepokeT(ip) ~= 0
            if metadata.TrialInfo_block2.nosepokeID(ip) == 5%add to nosepoke counts
                    for jk = 1:15001
                        if ip == 1
                            receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo_block2.trial_to_nosepokeT(ip)
                                receptacle1(jk,receptacle1_count) = NaN;
                            elseif temp / 1000 > metadata.TrialInfo_block2.nosepoke_length(ip)
                                receptacle1(jk:15001,receptacle1_count) = NaN;
                                break
                            else
                                receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    receptacle1_count = receptacle1_count + 1;

            end
            temp = -5000;
            nosepoke_number = nosepoke_number + 1;
        end
    end
            
        case 1 %R060
           for ik = 1:length(metadata.TrialInfo{1}.nosepokeT)
        if metadata.TrialInfo{1}.nosepokeT(ik) ~= 0
            if metadata.TrialInfo{1}.nosepokeID(ik) == 5 %add to nosepoke counts
                    for jk = 1:15001
                        if ik == 1
                            receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{1}.trial_to_nosepokeT(ik)
                                receptacle1(jk,receptacle1_count) = NaN;
                            elseif temp / 1000 > metadata.TrialInfo{1}.nosepoke_length(ik)
                                receptacle1(jk:15001,receptacle1_count) = NaN;
                                break
                            else
                                receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    receptacle1_count = receptacle1_count + 1;

            end
            temp = -5000;
            nosepoke_number = nosepoke_number + 1;
        end
    end
    
    for ip = 1:length(metadata.TrialInfo{2}.nosepokeT)
        if metadata.TrialInfo{2}.nosepokeT(ip) ~= 0
            if metadata.TrialInfo{2}.nosepokeID(ip) == 5%add to nosepoke counts
                    for jk = 1:15001
                        if ip == 1
                            receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                        else
                            if temp / -1000 > metadata.TrialInfo{2}.trial_to_nosepokeT(ip)
                                receptacle1(jk,receptacle1_count) = NaN;
                            elseif temp / 1000 > metadata.TrialInfo{2}.nosepoke_length(ip)
                                receptacle1(jk:15001,receptacle1_count) = NaN;
                                break
                            else
                                receptacle1(jk,receptacle1_count) = gau_sdf(dataPoint.Nosepokes(nosepoke_number) + temp);
                            end
                        end
                        temp = temp + 1;
                    end
                    receptacle1_count = receptacle1_count + 1;

            end
            temp = -5000;
            nosepoke_number = nosepoke_number + 1;
        end
    end
    end
    PETH.Receptacle.ALL.receptacle1 = receptacle1;
    
    
    %%
    rng('shuffle')
    num_trials = round(length(PETH.Receptacle.ALL.receptacle1(1,:)));
    temp_idx = randperm(num_trials);
    half_1 = temp_idx(1:num_trials/2);
    half_2 = temp_idx(num_trials/2+1:end);
    
    for iBin = 1:15001 %generate averaged responses
        if iBin <= length(PETH.Receptacle.ALL.receptacle1)
            PETHs.trials_receptacle1_1st_half{count}(iBin,:) = nanmean(PETH.Receptacle.ALL.receptacle1(iBin,half_1));
            PETHs.trials_receptacle1_2nd_half{count}(iBin,:) = nanmean(PETH.Receptacle.ALL.receptacle1(iBin,half_2));
        end
    end
    
    [maximum_receptacle1_1st(count),max_location_receptacle1_1st(count)] = max(PETHs.trials_receptacle1_1st_half{count}(cfg.start:cfg.end));
    [maximum_receptacle1_2nd(count),max_location_receptacle1_2nd(count)] = max(PETHs.trials_receptacle1_2nd_half{count}(cfg.start:cfg.end));
    
    count = count + 1;
    
    %         case 1
    %     end
end

[sortPETH.receptacle1_1st.sorted,sortPETH.receptacle1_1st.index_sort] = sort(max_location_receptacle1_1st);
[sortPETH.receptacle1_2nd.sorted,sortPETH.receptacle1_2nd.index_sort] = sort(max_location_receptacle1_2nd);

for jj = 1:length(mat_files_drift)
    sortedPETH.receptacle1_1st.unaltered(jj,:) = PETHs.trials_receptacle1_1st_half{sortPETH.receptacle1_1st.index_sort(jj)}(cfg.start:cfg.end);
    sortedPETH.receptacle1_2nd.receptacle1_1st(jj,:) = PETHs.trials_receptacle1_2nd_half{sortPETH.receptacle1_1st.index_sort(jj)}(cfg.start:cfg.end);
    sortedPETH.receptacle1_2nd.unaltered(jj,:) = PETHs.trials_receptacle1_2nd_half{sortPETH.receptacle1_2nd.index_sort(jj)}(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(mat_files_drift)
    mean_peth.receptacle1_1st(jj,1) = mean(sortedPETH.receptacle1_1st.unaltered(jj,:));
    std_peth.receptacle1_1st(jj,1) = std(sortedPETH.receptacle1_1st.unaltered(jj,:));
    sortedPETH.receptacle1_1st.zscore(jj,:) = (sortedPETH.receptacle1_1st.unaltered(jj,:) - mean_peth.receptacle1_1st(jj)) / std_peth.receptacle1_1st(jj);
    mean_peth.receptacle1_2nd(jj,1) = mean(sortedPETH.receptacle1_2nd.unaltered(jj,:));
    std_peth.receptacle1_2nd(jj,1) = std(sortedPETH.receptacle1_2nd.unaltered(jj,:));
    sortedPETH.receptacle1_2nd.zscore(jj,:) = (sortedPETH.receptacle1_2nd.unaltered(jj,:) - mean_peth.receptacle1_2nd(jj)) / std_peth.receptacle1_2nd(jj);
    mean_peth.receptacle1_2nd_2_1st(jj,1) = mean(sortedPETH.receptacle1_2nd.receptacle1_1st(jj,:));
    std_peth.receptacle1_2nd_2_1st(jj,1) = std(sortedPETH.receptacle1_2nd.receptacle1_1st(jj,:));
    sortedPETH.receptacle1_2nd_2_1st.zscore(jj,:) = (sortedPETH.receptacle1_2nd.receptacle1_1st(jj,:) - mean_peth.receptacle1_2nd_2_1st(jj)) / std_peth.receptacle1_2nd_2_1st(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','sortPETH','sortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:1,1:length(dir('*.mat')),sortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(mat_files_drift),sortedPETH.receptacle1_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(mat_files_drift),sortedPETH.receptacle1_2nd_2_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(mat_files_drift),sortedPETH.receptacle1_2nd.zscore);
colorbar; caxis([-3 4]);