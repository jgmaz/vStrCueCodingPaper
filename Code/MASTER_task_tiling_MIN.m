% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

cfg.start = 4001;
cfg.end = 8000;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    [minimum_light(kk),min_location_light(kk)] = min(PETH.Nosepoke.MEAN.nosepoke_light_PETH(cfg.start:cfg.end));
    [minimum_sound(kk),min_location_sound(kk)] = min(PETH.Nosepoke.MEAN.nosepoke_sound_PETH(cfg.start:cfg.end));
end

[MINsortPETH.light.sorted,MINsortPETH.light.index_sort] = sort(min_location_light);
[MINsortPETH.sound.sorted,MINsortPETH.sound.index_sort] = sort(min_location_sound);

for jj = 1:length(dir('*.mat'))
    load(mat_files(MINsortPETH.light.index_sort(jj)).name);
    MINsortedPETH.light.unaltered(jj,:) = PETH.Nosepoke.MEAN.nosepoke_light_PETH(cfg.start:cfg.end);
    MINsortedPETH.sound.light(jj,:) = PETH.Nosepoke.MEAN.nosepoke_sound_PETH(cfg.start:cfg.end);
    load(mat_files(MINsortPETH.sound.index_sort(jj)).name);
    MINsortedPETH.sound.unaltered(jj,:) = PETH.Nosepoke.MEAN.nosepoke_sound_PETH(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(dir('*.mat'))
    mean_peth.light(jj,1) = mean(MINsortedPETH.light.unaltered(jj,:));
    std_peth.light(jj,1) = std(MINsortedPETH.light.unaltered(jj,:));
    MINsortedPETH.light.zscore(jj,:) = (MINsortedPETH.light.unaltered(jj,:) - mean_peth.light(jj)) / std_peth.light(jj);
        mean_peth.sound(jj,1) = mean(MINsortedPETH.sound.unaltered(jj,:));
    std_peth.sound(jj,1) = std(MINsortedPETH.sound.unaltered(jj,:));
    MINsortedPETH.sound.zscore(jj,:) = (MINsortedPETH.sound.unaltered(jj,:) - mean_peth.sound(jj)) / std_peth.sound(jj);
        mean_peth.sound_2_light(jj,1) = mean(MINsortedPETH.sound.light(jj,:));
    std_peth.sound_2_light(jj,1) = std(MINsortedPETH.sound.light(jj,:));
    MINsortedPETH.sound_2_light.zscore(jj,:) = (MINsortedPETH.sound.light(jj,:) - mean_peth.sound_2_light(jj)) / std_peth.sound_2_light(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','MINsortPETH','MINsortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:1,1:length(dir('*.mat')),MINsortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.light.zscore); colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound_2_light.zscore);
 colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound.zscore);
colorbar; caxis([-3 4]);