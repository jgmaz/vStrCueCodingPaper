% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

cfg.start = 4001;
cfg.end = 6000;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    [maximum_light(kk),max_location_light(kk)] = max(PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end));
    [maximum_sound(kk),max_location_sound(kk)] = max(PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end));
end

[sortPETH.light.sorted,sortPETH.light.index_sort] = sort(max_location_light);
[sortPETH.sound.sorted,sortPETH.sound.index_sort] = sort(max_location_sound);

for jj = 1:length(dir('*.mat'))
    load(mat_files(sortPETH.light.index_sort(jj)).name);
    sortedPETH.light.unaltered(jj,:) = PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end);
    sortedPETH.sound.light(jj,:) = PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end);
    load(mat_files(sortPETH.sound.index_sort(jj)).name);
    sortedPETH.sound.unaltered(jj,:) = PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(dir('*.mat'))
    mean_peth.light(jj,1) = mean(sortedPETH.light.unaltered(jj,:));
    std_peth.light(jj,1) = std(sortedPETH.light.unaltered(jj,:));
    sortedPETH.light.zscore(jj,:) = (sortedPETH.light.unaltered(jj,:) - mean_peth.light(jj)) / std_peth.light(jj);
        mean_peth.sound(jj,1) = mean(sortedPETH.sound.unaltered(jj,:));
    std_peth.sound(jj,1) = std(sortedPETH.sound.unaltered(jj,:));
    sortedPETH.sound.zscore(jj,:) = (sortedPETH.sound.unaltered(jj,:) - mean_peth.sound(jj)) / std_peth.sound(jj);
        mean_peth.sound_2_light(jj,1) = mean(sortedPETH.sound.light(jj,:));
    std_peth.sound_2_light(jj,1) = std(sortedPETH.sound.light(jj,:));
    sortedPETH.sound_2_light.zscore(jj,:) = (sortedPETH.sound.light(jj,:) - mean_peth.sound_2_light(jj)) / std_peth.sound_2_light(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','sortPETH','sortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:1,1:length(dir('*.mat')),sortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.light.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound_2_light.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound.zscore);
colorbar; caxis([-3 4]);