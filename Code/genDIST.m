function TaskDist = genDIST(directory,destination,cfg,direction,feature,epoch)
% function TaskDist = genDIST(directory,destination)
%
%
% INPUTS:
%
% OUTPUTS:

cd(directory)

order_direction = {'max' 'min'};
cue_feature = {'identity' 'location' 'outcome'};
condition1 = {'light' 'arm1' 'reward'};
condition2 = {'sound' 'arm2' 'unreward'};
Epoch = {'cueonset' 'nosepoke'};

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    switch direction
        case 1
    [point1(kk),location1(kk)] = max(PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end));
    [point2(kk),location2(kk)] = max(PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end));
     case 2
    [point1(kk),location1(kk)] = min(PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end));
    [point2(kk),location2(kk)] = min(PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end));
    end
end

[sortPETH.(condition1{feature}).sorted,sortPETH.(condition1{feature}).index_sort] = sort(location1);
[sortPETH.(condition2{feature}).sorted,sortPETH.(condition2{feature}).index_sort] = sort(location2);

for jj = 1:length(dir('*.mat'))
    load(mat_files(sortPETH.(condition1{feature}).index_sort(jj)).name);
    sortedPETH.(condition1{feature}).unaltered(jj,:) = PETH.Trial.MEAN.trials_light_PETH(cfg.start:cfg.end);
    sortedPETH.(condition2{feature}).(condition1{feature})(jj,:) = PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end);
    load(mat_files(sortPETH.(condition2{feature}).index_sort(jj)).name);
    sortedPETH.(condition2{feature}).unaltered(jj,:) = PETH.Trial.MEAN.trials_sound_PETH(cfg.start:cfg.end);
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

save(cat(2,destination,'Distributed_coding_',Epoch{epoch},'_',cue_feature{feature},'_',order_direction{direction},'.mat'),'sortedPETH')

end