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
condition1 = {'light' 'arm1' 'rew'};
condition2 = {'sound' 'arm2' 'unrew'};
Epoch = {'Trial' 'Nosepoke'};
epoch_lower = {'trials' 'nosepoke'};

PETH_variable{1} = strcat(epoch_lower{epoch},'_',condition1{feature},'_PETH');
PETH_variable{2} = strcat(epoch_lower{epoch},'_',condition2{feature},'_PETH');

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    switch direction
        case 1
            [point1(kk),location1(kk)] = max(PETH.(Epoch{epoch}).MEAN.(PETH_variable{1})(cfg.start:cfg.end));
            [point2(kk),location2(kk)] = max(PETH.(Epoch{epoch}).MEAN.(PETH_variable{2})(cfg.start:cfg.end));
        case 2
            [point1(kk),location1(kk)] = min(PETH.(Epoch{epoch}).MEAN.(PETH_variable{1})(cfg.start:cfg.end));
            [point2(kk),location2(kk)] = min(PETH.(Epoch{epoch}).MEAN.(PETH_variable{2})(cfg.start:cfg.end));
    end
end

[sortPETH.(condition1{feature}).sorted,sortPETH.(condition1{feature}).index_sort] = sort(location1);
[sortPETH.(condition2{feature}).sorted,sortPETH.(condition2{feature}).index_sort] = sort(location2);

for jj = 1:length(dir('*.mat'))
    load(mat_files(sortPETH.(condition1{feature}).index_sort(jj)).name);
    sortedPETH.(condition1{feature}).unaltered(jj,:) = PETH.(Epoch{epoch}).MEAN.(PETH_variable{1})(cfg.start:cfg.end);
    sortedPETH.(condition2{feature}).(condition1{feature})(jj,:) = PETH.(Epoch{epoch}).MEAN.(PETH_variable{2})(cfg.start:cfg.end);
    load(mat_files(sortPETH.(condition2{feature}).index_sort(jj)).name);
    sortedPETH.(condition2{feature}).unaltered(jj,:) = PETH.(Epoch{epoch}).MEAN.(PETH_variable{2})(cfg.start:cfg.end);
    disp(jj)
end

for jj = 1:length(dir('*.mat'))
    mean_peth.(condition1{feature})(jj,1) = mean(sortedPETH.(condition1{feature}).unaltered(jj,:));
    std_peth.(condition1{feature})(jj,1) = std(sortedPETH.(condition1{feature}).unaltered(jj,:));
    sortedPETH.(condition1{feature}).zscore(jj,:) = (sortedPETH.(condition1{feature}).unaltered(jj,:) - mean_peth.(condition1{feature})(jj)) / std_peth.(condition1{feature})(jj);
    mean_peth.(condition2{feature})(jj,1) = mean(sortedPETH.(condition2{feature}).unaltered(jj,:));
    std_peth.(condition2{feature})(jj,1) = std(sortedPETH.(condition2{feature}).unaltered(jj,:));
    sortedPETH.(condition2{feature}).zscore(jj,:) = (sortedPETH.(condition2{feature}).unaltered(jj,:) - mean_peth.(condition2{feature})(jj)) / std_peth.(condition2{feature})(jj);
    mean_peth.TwoVsOne(jj,1) = mean(sortedPETH.(condition2{feature}).(condition1{feature})(jj,:));
    std_peth.TwoVsOne(jj,1) = std(sortedPETH.(condition2{feature}).(condition1{feature})(jj,:));
    sortedPETH.TwoVsOne.zscore(jj,:) = (sortedPETH.(condition2{feature}).(condition1{feature})(jj,:) - mean_peth.TwoVsOne(jj)) / std_peth.TwoVsOne(jj);
end

save(cat(2,destination,'Distributed_coding_',epoch_lower{epoch},'_',cue_feature{feature},'_',order_direction{direction},'.mat'),'sortedPETH')

end