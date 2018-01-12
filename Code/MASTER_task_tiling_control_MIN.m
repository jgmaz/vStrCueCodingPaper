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
    
    rng('shuffle')                
    num_trials = round(length(PETH.Trial.ALL.trials_rew_PETH(1,:)));
    temp_idx = randperm(num_trials);
    half_1 = temp_idx(1:num_trials/2);
    half_2 = temp_idx(num_trials/2+1:end);
    
    for iBin = 1:15001 %generate averaged responses
        if iBin <= length(PETH.Trial.ALL.trials_rew_PETH)
        PETHs.trials_rew_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_rew_PETH(iBin,half_1));
    PETHs.trials_rew_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_rew_PETH(iBin,half_2));
        end
    end
    
    [minimum_rew_1st(count),min_location_rew_1st(count)] = min(PETHs.trials_rew_1st_half{count}(cfg.start:cfg.end));
    [minimum_rew_2nd(count),min_location_rew_2nd(count)] = min(PETHs.trials_rew_2nd_half{count}(cfg.start:cfg.end));
    
    count = count + 1;
    
%         case 1
%     end
end

[MINsortPETH.rew_1st.sorted,MINsortPETH.rew_1st.index_sort] = sort(min_location_rew_1st);
[MINsortPETH.rew_2nd.sorted,MINsortPETH.rew_2nd.index_sort] = sort(min_location_rew_2nd);

for jj = 1:length(mat_files_drift)  
    MINsortedPETH.rew_1st.unaltered(jj,:) = PETHs.trials_rew_1st_half{MINsortPETH.rew_1st.index_sort(jj)}(cfg.start:cfg.end);
    MINsortedPETH.rew_2nd.rew_1st(jj,:) = PETHs.trials_rew_2nd_half{MINsortPETH.rew_1st.index_sort(jj)}(cfg.start:cfg.end);
    MINsortedPETH.rew_2nd.unaltered(jj,:) = PETHs.trials_rew_2nd_half{MINsortPETH.rew_2nd.index_sort(jj)}(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(mat_files_drift)
    mean_peth.rew_1st(jj,1) = mean(MINsortedPETH.rew_1st.unaltered(jj,:));
    std_peth.rew_1st(jj,1) = std(MINsortedPETH.rew_1st.unaltered(jj,:));
    MINsortedPETH.rew_1st.zscore(jj,:) = (MINsortedPETH.rew_1st.unaltered(jj,:) - mean_peth.rew_1st(jj)) / std_peth.rew_1st(jj);
        mean_peth.rew_2nd(jj,1) = mean(MINsortedPETH.rew_2nd.unaltered(jj,:));
    std_peth.rew_2nd(jj,1) = std(MINsortedPETH.rew_2nd.unaltered(jj,:));
    MINsortedPETH.rew_2nd.zscore(jj,:) = (MINsortedPETH.rew_2nd.unaltered(jj,:) - mean_peth.rew_2nd(jj)) / std_peth.rew_2nd(jj);
        mean_peth.rew_2nd_2_1st(jj,1) = mean(MINsortedPETH.rew_2nd.rew_1st(jj,:));
    std_peth.rew_2nd_2_1st(jj,1) = std(MINsortedPETH.rew_2nd.rew_1st(jj,:));
    MINsortedPETH.rew_2nd_2_1st.zscore(jj,:) = (MINsortedPETH.rew_2nd.rew_1st(jj,:) - mean_peth.rew_2nd_2_1st(jj)) / std_peth.rew_2nd_2_1st(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','MINsortPETH','MINsortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:1,1:length(dir('*.mat')),MINsortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.rew_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.rew_2nd_2_1st.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(mat_files_drift),MINsortedPETH.rew_2nd.zscore);
colorbar; caxis([-3 4]);