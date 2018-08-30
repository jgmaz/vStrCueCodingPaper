function DISTplot = plotDIST(directory,which_plot)
% function DISTplot = plotDIST(directory,which_plot)
%
%
% INPUTS:
%
% OUTPUTS:

cd(directory)

order_direction = {'max' 'min' 'max' 'min'};
cue_feature = {'identity' 'location' 'outcome'};
load_var = {'cueon' 'cueon' 'NP' 'NP'};
start_pt = [1 5 9];

if which_plot == 1 || which_plot == 2
    condition1 = {'light' 'photosensor1' 'rew'};
    condition2 = {'sound' 'photosensor2' 'unrew'};
else
    condition1 = {'light' 'receptacle1' 'rew'};
    condition2 = {'sound' 'receptacle2' 'unrew'};
end
figure('units','normalized','outerposition',[0 0 1 1]);
for iFeature = 1:length(cue_feature)
    
    switch iFeature
        case 1
            titles = {'Light block' 'Sound block' 'Sound vs. Light block'};
        case 2
            titles = {'Arm 1' 'Arm 2' 'Arm 2 vs. Arm 1'};
        case 3
            titles = {'Reward-available' 'Reward-unavailable' 'Reward-unavailable vs. -available'};
    end
    
    load(strcat(directory,'Distributed_coding_',load_var{which_plot},'_',cue_feature{iFeature},'_',order_direction{which_plot},'.mat'));
    
    Exclude = [];
    Exclude = find(isnan(sortedPETH.(condition1{iFeature}).zscore(:,1)));
    sortedPETH.(condition1{iFeature}).zscore(Exclude,:) = [];
    Exclude = [];
    Exclude = find(isnan(sortedPETH.(condition2{iFeature}).zscore(:,1)));
    sortedPETH.(condition2{iFeature}).zscore(Exclude,:) = [];
    Exclude = [];
    Exclude = find(isnan(sortedPETH.TwoVsOne.zscore(:,1)));
    sortedPETH.TwoVsOne.zscore(Exclude,:) = [];
    
    axis_length = length(sortedPETH.TwoVsOne.zscore)/1000 - 1;
    
    %%
    subplot(3,4,start_pt(iFeature))
    imagesc(-1:.001:axis_length,1:length(dir('*.mat')),sortedPETH.(condition1{iFeature}).zscore);
    hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
    caxis([-3 4]);
    title(titles{1});
    ylabel('Unit number');
    if iFeature == 3
        set(gca,'FontSize',18,'YTick',[])
        xlabel('Time from cue onset (s)');
    else
        set(gca,'FontSize',18,'YTick',[],'XTick',[]);
    end
    
    subplot(3,4,start_pt(iFeature)+1)
    imagesc(-1:.001:axis_length,1:length(dir('*.mat')),sortedPETH.(condition2{iFeature}).zscore);
    hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
    caxis([-3 4]);
    title(titles{2});
    if iFeature == 3
        set(gca,'FontSize',18,'YTick',[])
        xlabel('Time from cue onset (s)');
    else
        set(gca,'FontSize',18,'YTick',[],'XTick',[]);
    end
    
    subplot(3,4,start_pt(iFeature)+2)
    imagesc(-1:.001:1,axis_length:length(dir('*.mat')),sortedPETH.TwoVsOne.zscore);
    hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
    caxis([-3 4]);
    title(titles{3});
    if iFeature == 3
        set(gca,'FontSize',18,'YTick',[])
        xlabel('Time from cue onset (s)');
    else
        set(gca,'FontSize',18,'YTick',[],'XTick',[]);
    end
    
end