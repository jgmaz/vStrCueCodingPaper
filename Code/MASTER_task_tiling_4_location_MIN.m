% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

cfg.start = 4001;
cfg.end = 8000;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    [minimum_arm1(kk),min_location_arm1(kk)] = min(PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end));
    [minimum_arm2(kk),min_location_arm2(kk)] = min(PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end));
    [minimum_arm3(kk),min_location_arm3(kk)] = min(PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end));
    [minimum_arm4(kk),min_location_arm4(kk)] = min(PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end));
end

[MINsortPETH.arm1.sorted,MINsortPETH.arm1.index_sort] = sort(min_location_arm1);
[MINsortPETH.arm2.sorted,MINsortPETH.arm2.index_sort] = sort(min_location_arm2);
[MINsortPETH.arm3.sorted,MINsortPETH.arm3.index_sort] = sort(min_location_arm3);
[MINsortPETH.arm4.sorted,MINsortPETH.arm4.index_sort] = sort(min_location_arm4);

for jj = 1:length(dir('*.mat'))
    load(mat_files(MINsortPETH.arm1.index_sort(jj)).name);
    MINsortedPETH.arm1.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end);
    MINsortedPETH.arm2.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end);
    MINsortedPETH.arm3.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    MINsortedPETH.arm4.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    load(mat_files(MINsortPETH.arm2.index_sort(jj)).name);
    MINsortedPETH.arm2.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end);
    MINsortedPETH.arm3.arm2(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    load(mat_files(MINsortPETH.arm3.index_sort(jj)).name);
    MINsortedPETH.arm3.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    MINsortedPETH.arm4.arm3(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    load(mat_files(MINsortPETH.arm4.index_sort(jj)).name);
    MINsortedPETH.arm4.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    MINsortedPETH.arm1.arm4(jj,:) = PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(dir('*.mat'))
    mean_peth.arm1(jj,1) = mean(MINsortedPETH.arm1.unaltered(jj,:));
    std_peth.arm1(jj,1) = std(MINsortedPETH.arm1.unaltered(jj,:));
    MINsortedPETH.arm1.zscore(jj,:) = (MINsortedPETH.arm1.unaltered(jj,:) - mean_peth.arm1(jj)) / std_peth.arm1(jj);
    mean_peth.arm1_2_arm4(jj,1) = mean(MINsortedPETH.arm1.arm4(jj,:));
    std_peth.arm1_2_arm4(jj,1) = std(MINsortedPETH.arm1.arm4(jj,:));
    MINsortedPETH.arm1_2_arm4.zscore(jj,:) = (MINsortedPETH.arm1.arm4(jj,:) - mean_peth.arm1_2_arm4(jj)) / std_peth.arm1_2_arm4(jj);
    mean_peth.arm2(jj,1) = mean(MINsortedPETH.arm2.unaltered(jj,:));
    std_peth.arm2(jj,1) = std(MINsortedPETH.arm2.unaltered(jj,:));
    MINsortedPETH.arm2.zscore(jj,:) = (MINsortedPETH.arm2.unaltered(jj,:) - mean_peth.arm2(jj)) / std_peth.arm2(jj);
    mean_peth.arm2_2_arm1(jj,1) = mean(MINsortedPETH.arm2.arm1(jj,:));
    std_peth.arm2_2_arm1(jj,1) = std(MINsortedPETH.arm2.arm1(jj,:));
    MINsortedPETH.arm2_2_arm1.zscore(jj,:) = (MINsortedPETH.arm2.arm1(jj,:) - mean_peth.arm2_2_arm1(jj)) / std_peth.arm2_2_arm1(jj);
    mean_peth.arm3(jj,1) = mean(MINsortedPETH.arm3.unaltered(jj,:));
    std_peth.arm3(jj,1) = std(MINsortedPETH.arm3.unaltered(jj,:));
    MINsortedPETH.arm3.zscore(jj,:) = (MINsortedPETH.arm3.unaltered(jj,:) - mean_peth.arm3(jj)) / std_peth.arm3(jj);
    mean_peth.arm3_2_arm1(jj,1) = mean(MINsortedPETH.arm3.arm1(jj,:));
    std_peth.arm3_2_arm1(jj,1) = std(MINsortedPETH.arm3.arm1(jj,:));
    MINsortedPETH.arm3_2_arm1.zscore(jj,:) = (MINsortedPETH.arm3.arm1(jj,:) - mean_peth.arm3_2_arm1(jj)) / std_peth.arm3_2_arm1(jj);
    mean_peth.arm3_2_arm2(jj,1) = mean(MINsortedPETH.arm3.arm2(jj,:));
    std_peth.arm3_2_arm2(jj,1) = std(MINsortedPETH.arm3.arm2(jj,:));
    MINsortedPETH.arm3_2_arm2.zscore(jj,:) = (MINsortedPETH.arm3.arm2(jj,:) - mean_peth.arm3_2_arm2(jj)) / std_peth.arm3_2_arm2(jj);
    mean_peth.arm4(jj,1) = mean(MINsortedPETH.arm4.unaltered(jj,:));
    std_peth.arm4(jj,1) = std(MINsortedPETH.arm4.unaltered(jj,:));
    MINsortedPETH.arm4.zscore(jj,:) = (MINsortedPETH.arm4.unaltered(jj,:) - mean_peth.arm4(jj)) / std_peth.arm4(jj);
    mean_peth.arm4_2_arm1(jj,1) = mean(MINsortedPETH.arm4.arm1(jj,:));
    std_peth.arm4_2_arm1(jj,1) = std(MINsortedPETH.arm4.arm1(jj,:));
    MINsortedPETH.arm4_2_arm1.zscore(jj,:) = (MINsortedPETH.arm4.arm1(jj,:) - mean_peth.arm4_2_arm1(jj)) / std_peth.arm4_2_arm1(jj);
    mean_peth.arm4_2_arm3(jj,1) = mean(MINsortedPETH.arm4.arm3(jj,:));
    std_peth.arm4_2_arm3(jj,1) = std(MINsortedPETH.arm4.arm3(jj,:));
    MINsortedPETH.arm4_2_arm3.zscore(jj,:) = (MINsortedPETH.arm4.arm3(jj,:) - mean_peth.arm4_2_arm3(jj)) / std_peth.arm4_2_arm3(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','MINsortPETH','MINsortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm1.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2_2_arm1.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm1_2_arm4.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm3_2_arm2.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm4.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm4_2_arm3.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm3.zscore);
colorbar; caxis([-3 4]);