% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

cfg.start = 4001;
cfg.end = 6000;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    [maximum_arm1(kk),max_location_arm1(kk)] = max(PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end));
    [maximum_arm2(kk),max_location_arm2(kk)] = max(PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end));
    [maximum_arm3(kk),max_location_arm3(kk)] = max(PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end));
    [maximum_arm4(kk),max_location_arm4(kk)] = max(PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end));
end

[sortPETH.arm1.sorted,sortPETH.arm1.index_sort] = sort(max_location_arm1);
[sortPETH.arm2.sorted,sortPETH.arm2.index_sort] = sort(max_location_arm2);
[sortPETH.arm3.sorted,sortPETH.arm3.index_sort] = sort(max_location_arm3);
[sortPETH.arm4.sorted,sortPETH.arm4.index_sort] = sort(max_location_arm4);

for jj = 1:length(dir('*.mat'))
    load(mat_files(sortPETH.arm1.index_sort(jj)).name);
    sortedPETH.arm1.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end);
    sortedPETH.arm2.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end);
    sortedPETH.arm3.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    sortedPETH.arm4.arm1(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    load(mat_files(sortPETH.arm2.index_sort(jj)).name);
    sortedPETH.arm2.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle2(cfg.start:cfg.end);
    sortedPETH.arm3.arm2(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    load(mat_files(sortPETH.arm3.index_sort(jj)).name);
    sortedPETH.arm3.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle3(cfg.start:cfg.end);
    sortedPETH.arm4.arm3(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    load(mat_files(sortPETH.arm4.index_sort(jj)).name);
    sortedPETH.arm4.unaltered(jj,:) = PETH.Receptacle.MEAN.receptacle4(cfg.start:cfg.end);
    sortedPETH.arm1.arm4(jj,:) = PETH.Receptacle.MEAN.receptacle1(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(dir('*.mat'))
    mean_peth.arm1(jj,1) = mean(sortedPETH.arm1.unaltered(jj,:));
    std_peth.arm1(jj,1) = std(sortedPETH.arm1.unaltered(jj,:));
    sortedPETH.arm1.zscore(jj,:) = (sortedPETH.arm1.unaltered(jj,:) - mean_peth.arm1(jj)) / std_peth.arm1(jj);
    mean_peth.arm1_2_arm4(jj,1) = mean(sortedPETH.arm1.arm4(jj,:));
    std_peth.arm1_2_arm4(jj,1) = std(sortedPETH.arm1.arm4(jj,:));
    sortedPETH.arm1_2_arm4.zscore(jj,:) = (sortedPETH.arm1.arm4(jj,:) - mean_peth.arm1_2_arm4(jj)) / std_peth.arm1_2_arm4(jj);
    mean_peth.arm2(jj,1) = mean(sortedPETH.arm2.unaltered(jj,:));
    std_peth.arm2(jj,1) = std(sortedPETH.arm2.unaltered(jj,:));
    sortedPETH.arm2.zscore(jj,:) = (sortedPETH.arm2.unaltered(jj,:) - mean_peth.arm2(jj)) / std_peth.arm2(jj);
    mean_peth.arm2_2_arm1(jj,1) = mean(sortedPETH.arm2.arm1(jj,:));
    std_peth.arm2_2_arm1(jj,1) = std(sortedPETH.arm2.arm1(jj,:));
    sortedPETH.arm2_2_arm1.zscore(jj,:) = (sortedPETH.arm2.arm1(jj,:) - mean_peth.arm2_2_arm1(jj)) / std_peth.arm2_2_arm1(jj);
    mean_peth.arm3(jj,1) = mean(sortedPETH.arm3.unaltered(jj,:));
    std_peth.arm3(jj,1) = std(sortedPETH.arm3.unaltered(jj,:));
    sortedPETH.arm3.zscore(jj,:) = (sortedPETH.arm3.unaltered(jj,:) - mean_peth.arm3(jj)) / std_peth.arm3(jj);
    mean_peth.arm3_2_arm1(jj,1) = mean(sortedPETH.arm3.arm1(jj,:));
    std_peth.arm3_2_arm1(jj,1) = std(sortedPETH.arm3.arm1(jj,:));
    sortedPETH.arm3_2_arm1.zscore(jj,:) = (sortedPETH.arm3.arm1(jj,:) - mean_peth.arm3_2_arm1(jj)) / std_peth.arm3_2_arm1(jj);
    mean_peth.arm3_2_arm2(jj,1) = mean(sortedPETH.arm3.arm2(jj,:));
    std_peth.arm3_2_arm2(jj,1) = std(sortedPETH.arm3.arm2(jj,:));
    sortedPETH.arm3_2_arm2.zscore(jj,:) = (sortedPETH.arm3.arm2(jj,:) - mean_peth.arm3_2_arm2(jj)) / std_peth.arm3_2_arm2(jj);
    mean_peth.arm4(jj,1) = mean(sortedPETH.arm4.unaltered(jj,:));
    std_peth.arm4(jj,1) = std(sortedPETH.arm4.unaltered(jj,:));
    sortedPETH.arm4.zscore(jj,:) = (sortedPETH.arm4.unaltered(jj,:) - mean_peth.arm4(jj)) / std_peth.arm4(jj);
    mean_peth.arm4_2_arm1(jj,1) = mean(sortedPETH.arm4.arm1(jj,:));
    std_peth.arm4_2_arm1(jj,1) = std(sortedPETH.arm4.arm1(jj,:));
    sortedPETH.arm4_2_arm1.zscore(jj,:) = (sortedPETH.arm4.arm1(jj,:) - mean_peth.arm4_2_arm1(jj)) / std_peth.arm4_2_arm1(jj);
    mean_peth.arm4_2_arm3(jj,1) = mean(sortedPETH.arm4.arm3(jj,:));
    std_peth.arm4_2_arm3(jj,1) = std(sortedPETH.arm4.arm3(jj,:));
    sortedPETH.arm4_2_arm3.zscore(jj,:) = (sortedPETH.arm4.arm3(jj,:) - mean_peth.arm4_2_arm3(jj)) / std_peth.arm4_2_arm3(jj);
end

% save(cat(2,'E:\Jimmie\Jimmie\Analysis\','sortPETH','sortedPETH'));
%% plot
% figure
% colormap('default');   % set colormap
% imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.unaltered);        % draw image and scale colormap to values range
% colorbar;
figure;
subtightplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm1.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2_2_arm1.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm1_2_arm4.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm3_2_arm2.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm4.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm4_2_arm3.zscore);
colorbar; caxis([-3 4]);
subtightplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm3.zscore);
colorbar; caxis([-3 4]);