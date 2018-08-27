%% rewarded vs unrewarded trials for both blocks
mat_files = dir('Events*');
for iFind = 1:length(mat_files)
    
    rew_trial_num_block1 = 1;
    unrew_trial_num_block1 = 1;
    rew_trial_num_block2 = 1;
    unrew_trial_num_block2 = 1;
    
    app_rew_trials_block1 = [];
    app_unrew_trials_block1 = [];
    app_rew_trials_block2 = [];
    app_unrew_trials_block2 = [];
    
    length_rew_trials_block1 = [];
    length_unrew_trials_block1 = [];
    length_rew_trials_block2 = [];
    length_unrew_trials_block2 = [];
    
    avg_rew_trials_block1 = [];
    SEM_rew_trials_block1 = [];
    avg_unrew_trials_block1 = [];
    SEM_unrew_trials_block1 = [];
    avg_rew_trials_block2 = [];
    SEM_rew_trials_block2 = [];
    avg_unrew_trials_block2 = [];
    SEM_unrew_trials_block2 = [];
    
    for ik = 1:length(metadata.TrialInfo{1,1}.trialT)
        switch metadata.TrialInfo{1,1}.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                app_rew_trials_block1(rew_trial_num_block1) = metadata.TrialInfo{1,1}.summary(ik,3);
                if metadata.TrialInfo{1,1}.summary(ik,15) < 5
                    length_rew_trials_block1(rew_trial_num_block1) = metadata.TrialInfo{1,1}.summary(ik,15);
                else
                    length_rew_trials_block1(rew_trial_num_block1) = NaN;
                end
                rew_trial_num_block1 = rew_trial_num_block1 + 1;
            case 0
                app_unrew_trials_block1(unrew_trial_num_block1) = metadata.TrialInfo{1,1}.summary(ik,3);
                if metadata.TrialInfo{1,1}.summary(ik,15) < 5
                    length_unrew_trials_block1(unrew_trial_num_block1) = metadata.TrialInfo{1,1}.summary(ik,15);
                else
                    length_unrew_trials_block1(unrew_trial_num_block1) = NaN;
                end
                unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
        end
    end
    
    for ip = 1:length(metadata.TrialInfo{1,2}.trialT)
        switch metadata.TrialInfo{1,2}.rewarded(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                app_rew_trials_block2(rew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,3);
                if metadata.TrialInfo{1,2}.summary(ip,15) < 5
                    length_rew_trials_block2(rew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,15);
                else
                    length_rew_trials_block2(rew_trial_num_block2) = NaN;
                end
                rew_trial_num_block2 = rew_trial_num_block2 + 1;
            case 0
                app_unrew_trials_block2(unrew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,3);
                if metadata.TrialInfo{1,2}.summary(ip,15) < 5
                    length_unrew_trials_block2(unrew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,15);
                else
                    length_unrew_trials_block2(unrew_trial_num_block2) = NaN;
                end
                unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
        end
    end
    %% summary
    
    count_rew_trials_block1 = length(app_rew_trials_block1);
    prop_app_rew_trials_block1 = sum(app_rew_trials_block1)/count_rew_trials_block1;
    avg_length_rew_trials_block1 = nanmean(length_rew_trials_block1);
    SEM_length_rew_trials_block1 = nanstd(length_rew_trials_block1/sqrt(numel(length_rew_trials_block1)-sum(isnan(length_rew_trials_block1))));
    
    count_unrew_trials_block1 = length(app_unrew_trials_block1);
    prop_app_unrew_trials_block1 = sum(app_unrew_trials_block1)/count_unrew_trials_block1;
    avg_length_unrew_trials_block1 = nanmean(length_unrew_trials_block1);
    SEM_length_unrew_trials_block1 = nanstd(length_unrew_trials_block1/sqrt(numel(length_unrew_trials_block1)-sum(isnan(length_unrew_trials_block1))));
    
    count_rew_trials_block2 = length(app_rew_trials_block2);
    prop_app_rew_trials_block2 = sum(app_rew_trials_block2)/count_rew_trials_block2;
    avg_length_rew_trials_block2 = nanmean(length_rew_trials_block2);
    SEM_length_rew_trials_block2 = nanstd(length_rew_trials_block2/sqrt(numel(length_rew_trials_block2)-sum(isnan(length_rew_trials_block2))));
    
    count_unrew_trials_block2 = length(app_unrew_trials_block2);
    prop_app_unrew_trials_block2 = sum(app_unrew_trials_block2)/count_unrew_trials_block2;
    avg_length_unrew_trials_block2 = nanmean(length_unrew_trials_block2);
    SEM_length_unrew_trials_block2 = nanstd(length_unrew_trials_block2/sqrt(numel(length_unrew_trials_block2)-sum(isnan(length_unrew_trials_block2))));
    
    %% output
    switch sesh.block_order
        case 1
            BEHAV.Length.MEAN.rew_trials_light = avg_length_rew_trials_block1;
            BEHAV.Length.MEAN.unrew_trials_light = avg_length_unrew_trials_block1;
            BEHAV.Length.MEAN.rew_trials_sound = avg_length_rew_trials_block2;
            BEHAV.Length.MEAN.unrew_trials_sound = avg_length_unrew_trials_block2;
            
            BEHAV.Length.SEM.rew_trials_light = SEM_length_rew_trials_block1;
            BEHAV.Length.SEM.unrew_trials_light = SEM_length_unrew_trials_block1;
            BEHAV.Length.SEM.rew_trials_sound = SEM_length_rew_trials_block2;
            BEHAV.Length.SEM.unrew_trials_sound = SEM_length_unrew_trials_block2;
            
            BEHAV.Length.ALL.rew_trials_light = length_rew_trials_block1;
            BEHAV.Length.ALL.unrew_trials_light = length_unrew_trials_block1;
            BEHAV.Length.ALL.rew_trials_sound = length_rew_trials_block2;
            BEHAV.Length.ALL.unrew_trials_sound = length_unrew_trials_block2;
            
            BEHAV.Approach.COUNT.rew_trials_light = count_rew_trials_block1;
            BEHAV.Approach.COUNT.unrew_trials_light = count_unrew_trials_block1;
            BEHAV.Approach.COUNT.rew_trials_sound = count_rew_trials_block2;
            BEHAV.Approach.COUNT.unrew_trials_sound = count_unrew_trials_block2;
            
            BEHAV.Approach.PROP.rew_trials_light = prop_app_rew_trials_block1;
            BEHAV.Approach.PROP.unrew_trials_light = prop_app_unrew_trials_block1;
            BEHAV.Approach.PROP.rew_trials_sound = prop_app_rew_trials_block2;
            BEHAV.Approach.PROP.unrew_trials_sound = prop_app_unrew_trials_block2;
            
            BEHAV.Approach.ALL.rew_trials_light = app_rew_trials_block1;
            BEHAV.Approach.ALL.unrew_trials_light = app_unrew_trials_block1;
            BEHAV.Approach.ALL.rew_trials_sound = app_rew_trials_block2;
            BEHAV.Approach.ALL.unrew_trials_sound = app_unrew_trials_block2;
            
        case 2
            BEHAV.Length.MEAN.rew_trials_light = avg_length_rew_trials_block2;
            BEHAV.Length.MEAN.unrew_trials_light = avg_length_unrew_trials_block2;
            BEHAV.Length.MEAN.rew_trials_sound = avg_length_rew_trials_block1;
            BEHAV.Length.MEAN.unrew_trials_sound = avg_length_unrew_trials_block1;
            
            BEHAV.Length.SEM.rew_trials_light = SEM_length_rew_trials_block2;
            BEHAV.Length.SEM.unrew_trials_light = SEM_length_unrew_trials_block2;
            BEHAV.Length.SEM.rew_trials_sound = SEM_length_rew_trials_block1;
            BEHAV.Length.SEM.unrew_trials_sound = SEM_length_unrew_trials_block1;
            
            BEHAV.Length.ALL.rew_trials_light = length_rew_trials_block2;
            BEHAV.Length.ALL.unrew_trials_light = length_unrew_trials_block2;
            BEHAV.Length.ALL.rew_trials_sound = length_rew_trials_block1;
            BEHAV.Length.ALL.unrew_trials_sound = length_unrew_trials_block1;
            
            BEHAV.Approach.COUNT.rew_trials_light = count_rew_trials_block2;
            BEHAV.Approach.COUNT.unrew_trials_light = count_unrew_trials_block2;
            BEHAV.Approach.COUNT.rew_trials_sound = count_rew_trials_block1;
            BEHAV.Approach.COUNT.unrew_trials_sound = count_unrew_trials_block1;
            
            BEHAV.Approach.PROP.rew_trials_light = prop_app_rew_trials_block2;
            BEHAV.Approach.PROP.unrew_trials_light = prop_app_unrew_trials_block2;
            BEHAV.Approach.PROP.rew_trials_sound = prop_app_rew_trials_block1;
            BEHAV.Approach.PROP.unrew_trials_sound = prop_app_unrew_trials_block1;
            
            BEHAV.Approach.ALL.rew_trials_light = app_rew_trials_block2;
            BEHAV.Approach.ALL.unrew_trials_light = app_unrew_trials_block2;
            BEHAV.Approach.ALL.rew_trials_sound = app_rew_trials_block1;
            BEHAV.Approach.ALL.unrew_trials_sound = app_unrew_trials_block1;
    end
    
    %%
    switch which_constraints
        case 1
            save(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat\BEHAV\BEHAV\',sesh.session_id),'BEHAV');
        case 2
            save(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat\BEHAV\BEHAV2\',sesh.session_id),'BEHAV');
        case 3
            save(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat\BEHAV\BEHAV3\',sesh.session_id),'BEHAV');
    end
    clearvars -except metadata sesh
    clear
    %% group
    % cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R053\Mat\BEHAV\BEHAV'));
    
    counter = 1;
    mat_files = dir('*.mat');
    for kk = 1:length(dir('*.mat'))
        disp(num2str(kk));
        load(mat_files(kk).name);
        mat_overview.fname{kk} = mat_files(kk).name;
        
        BEHAV_summary.APP.rew_trials_light(kk) =  BEHAV.Approach.PROP.rew_trials_light;
        BEHAV_summary.APP.unrew_trials_light(kk) =  BEHAV.Approach.PROP.unrew_trials_light;
        BEHAV_summary.APP.rew_trials_sound(kk) =  BEHAV.Approach.PROP.rew_trials_sound;
        BEHAV_summary.APP.unrew_trials_sound(kk) =  BEHAV.Approach.PROP.unrew_trials_sound;
        BEHAV_summary.Length.rew_trials_light(kk) =  BEHAV.Length.MEAN.rew_trials_light;
        BEHAV_summary.Length.unrew_trials_light(kk) =  BEHAV.Length.MEAN.unrew_trials_light;
        BEHAV_summary.Length.rew_trials_sound(kk) =  BEHAV.Length.MEAN.rew_trials_sound;
        BEHAV_summary.Length.unrew_trials_sound(kk) =  BEHAV.Length.MEAN.unrew_trials_sound;
    end
    
    %% group summary
    
    BEHAV_summary.APP.MEAN.rew_trials_light = mean(BEHAV_summary.APP.rew_trials_light);
    BEHAV_summary.APP.MEAN.unrew_trials_light = mean(BEHAV_summary.APP.unrew_trials_light);
    BEHAV_summary.APP.MEAN.rew_trials_sound = mean(BEHAV_summary.APP.rew_trials_sound);
    BEHAV_summary.APP.MEAN.unrew_trials_sound = mean(BEHAV_summary.APP.unrew_trials_sound);
    
    BEHAV_summary.Length.MEAN.rew_trials_light = mean(BEHAV_summary.Length.rew_trials_light);
    BEHAV_summary.Length.MEAN.unrew_trials_light = mean(BEHAV_summary.Length.unrew_trials_light);
    BEHAV_summary.Length.MEAN.rew_trials_sound = mean(BEHAV_summary.Length.rew_trials_sound);
    BEHAV_summary.Length.MEAN.unrew_trials_sound = mean(BEHAV_summary.Length.unrew_trials_sound);
    
    BEHAV_summary.APP.SEM.rew_trials_light = std(BEHAV_summary.APP.rew_trials_light)/sqrt(numel(BEHAV_summary.APP.rew_trials_light));
    BEHAV_summary.APP.SEM.unrew_trials_light = std(BEHAV_summary.APP.unrew_trials_light)/sqrt(numel(BEHAV_summary.APP.unrew_trials_light));
    BEHAV_summary.APP.SEM.rew_trials_sound = std(BEHAV_summary.APP.rew_trials_sound)/sqrt(numel(BEHAV_summary.APP.rew_trials_sound));
    BEHAV_summary.APP.SEM.unrew_trials_sound = std(BEHAV_summary.APP.unrew_trials_sound)/sqrt(numel(BEHAV_summary.APP.unrew_trials_sound));
    
    BEHAV_summary.Length.SEM.rew_trials_light = std(BEHAV_summary.Length.rew_trials_light)/sqrt(numel(BEHAV_summary.Length.rew_trials_light));
    BEHAV_summary.Length.SEM.unrew_trials_light = std(BEHAV_summary.Length.unrew_trials_light)/sqrt(numel(BEHAV_summary.Length.unrew_trials_light));
    BEHAV_summary.Length.SEM.rew_trials_sound = std(BEHAV_summary.Length.rew_trials_sound)/sqrt(numel(BEHAV_summary.Length.rew_trials_sound));
    BEHAV_summary.Length.SEM.unrew_trials_sound = std(BEHAV_summary.Length.unrew_trials_sound)/sqrt(numel(BEHAV_summary.Length.unrew_trials_sound));
    
    %%
    BEHAV_summary.APP.ALL(:,1) = BEHAV_summary.APP.rew_trials_light;
    BEHAV_summary.APP.ALL(:,2) = BEHAV_summary.APP.unrew_trials_light;
    BEHAV_summary.APP.ALL(:,3) = BEHAV_summary.APP.rew_trials_sound;
    BEHAV_summary.APP.ALL(:,4) = BEHAV_summary.APP.unrew_trials_sound;
    
    BEHAV_summary.Length.ALL(:,1) = BEHAV_summary.Length.rew_trials_light;
    BEHAV_summary.Length.ALL(:,2) = BEHAV_summary.Length.unrew_trials_light;
    BEHAV_summary.Length.ALL(:,3) = BEHAV_summary.Length.rew_trials_sound;
    BEHAV_summary.Length.ALL(:,4) = BEHAV_summary.Length.unrew_trials_sound;
    
    %%
    
    [p_APP,tbl_APP,stats_APP] = anova1(BEHAV_summary.APP.ALL);
    [p_Length,tbl_Length,stats_Length] = anova1(BEHAV_summary.Length.ALL);
    
    
    [h_rvr p_rvr] = ttest(BEHAV_summary.APP.rew_trials_sound,BEHAV_summary.APP.rew_trials_light);
    [h_uvu p_uvu] = ttest(BEHAV_summary.APP.unrew_trials_sound,BEHAV_summary.APP.unrew_trials_light);
    [h_rvu1 p_rvu1] = ttest(BEHAV_summary.APP.rew_trials_light,BEHAV_summary.APP.unrew_trials_light);
    [h_rvu2 p_rvu2] = ttest(BEHAV_summary.APP.rew_trials_sound,BEHAV_summary.APP.unrew_trials_sound);
    
    [lh_rvr lp_rvr] = ttest(BEHAV_summary.Length.rew_trials_sound,BEHAV_summary.Length.rew_trials_light);
    [lh_uvu lp_uvu] = ttest(BEHAV_summary.Length.unrew_trials_sound,BEHAV_summary.Length.unrew_trials_light);
    [lh_rvu1 lp_rvu1] = ttest(BEHAV_summary.Length.rew_trials_light,BEHAV_summary.Length.unrew_trials_light);
    [lh_rvu2 lp_rvu2] = ttest(BEHAV_summary.Length.rew_trials_sound,BEHAV_summary.Length.unrew_trials_sound);
    
    %% all rats together
    BEHAV_summary.APP.MEAN.rew_trials_light = cat(2,R053_BEHAV_summary.APP.MEAN.rew_trials_light, ...
        R056_BEHAV_summary.APP.MEAN.rew_trials_light,R057_BEHAV_summary.APP.MEAN.rew_trials_light,R060_BEHAV_summary.APP.MEAN.rew_trials_light);
    BEHAV_summary.APP.MEAN.unrew_trials_light = cat(2,R053_BEHAV_summary.APP.MEAN.unrew_trials_light, ...
        R056_BEHAV_summary.APP.MEAN.unrew_trials_light,R057_BEHAV_summary.APP.MEAN.unrew_trials_light,R060_BEHAV_summary.APP.MEAN.unrew_trials_light);
    BEHAV_summary.APP.MEAN.rew_trials_sound = cat(2,R053_BEHAV_summary.APP.MEAN.rew_trials_sound, ...
        R056_BEHAV_summary.APP.MEAN.rew_trials_sound,R057_BEHAV_summary.APP.MEAN.rew_trials_sound,R060_BEHAV_summary.APP.MEAN.rew_trials_sound);
    BEHAV_summary.APP.MEAN.unrew_trials_sound = cat(2,R053_BEHAV_summary.APP.MEAN.unrew_trials_sound, ...
        R056_BEHAV_summary.APP.MEAN.unrew_trials_sound,R057_BEHAV_summary.APP.MEAN.unrew_trials_sound,R060_BEHAV_summary.APP.MEAN.unrew_trials_sound);
    
    BEHAV_summary.Length.MEAN.rew_trials_light = cat(2,R053_BEHAV_summary.Length.MEAN.rew_trials_light, ...
        R056_BEHAV_summary.Length.MEAN.rew_trials_light,R057_BEHAV_summary.Length.MEAN.rew_trials_light,R060_BEHAV_summary.Length.MEAN.rew_trials_light);
    BEHAV_summary.Length.MEAN.unrew_trials_light = cat(2,R053_BEHAV_summary.Length.MEAN.unrew_trials_light, ...
        R056_BEHAV_summary.Length.MEAN.unrew_trials_light,R057_BEHAV_summary.Length.MEAN.unrew_trials_light,R060_BEHAV_summary.Length.MEAN.unrew_trials_light);
    BEHAV_summary.Length.MEAN.rew_trials_sound = cat(2,R053_BEHAV_summary.Length.MEAN.rew_trials_sound, ...
        R056_BEHAV_summary.Length.MEAN.rew_trials_sound,R057_BEHAV_summary.Length.MEAN.rew_trials_sound,R060_BEHAV_summary.Length.MEAN.rew_trials_sound);
    BEHAV_summary.Length.MEAN.unrew_trials_sound = cat(2,R053_BEHAV_summary.Length.MEAN.unrew_trials_sound, ...
        R056_BEHAV_summary.Length.MEAN.unrew_trials_sound,R057_BEHAV_summary.Length.MEAN.unrew_trials_sound,R060_BEHAV_summary.Length.MEAN.unrew_trials_sound);
    
    %%
    mean_rew_light = mean(BEHAV_summary.APP.MEAN.rew_trials_light);
    mean_unrew_light = mean(BEHAV_summary.APP.MEAN.unrew_trials_light);
    mean_rew_sound = mean(BEHAV_summary.APP.MEAN.rew_trials_sound);
    mean_unrew_sound = mean(BEHAV_summary.APP.MEAN.unrew_trials_sound);
    
    mean_rew_light2 = mean(BEHAV_summary.Length.MEAN.rew_trials_light);
    mean_unrew_light2 = mean(BEHAV_summary.Length.MEAN.unrew_trials_light);
    mean_rew_sound2 = mean(BEHAV_summary.Length.MEAN.rew_trials_sound);
    mean_unrew_sound2 = mean(BEHAV_summary.Length.MEAN.unrew_trials_sound);
    %%
    group = [1 2 3 4];
    figure
    subplot(2,2,[1 3])
    gscatter([1 1 1 1],BEHAV_summary.APP.MEAN.rew_trials_light,group,'r','xo+*')
    hold on;
    gscatter([2 2 2 2],BEHAV_summary.APP.MEAN.unrew_trials_light,group,'g','xo+*')
    gscatter([3 3 3 3],BEHAV_summary.APP.MEAN.rew_trials_sound,group,'b','xo+*')
    gscatter([4 4 4 4],BEHAV_summary.APP.MEAN.unrew_trials_sound,group,'c','xo+*')
    plot([.85 1.15],[mean_rew_light mean_rew_light],'k')
    plot([1.85 2.15],[mean_unrew_light mean_unrew_light],'k')
    plot([2.85 3.15],[mean_rew_sound mean_rew_sound],'k')
    plot([3.85 4.15],[mean_unrew_sound mean_unrew_sound],'k')
    
    xlim([0 5]); ylim([0 1]); title('Proportion Approached');
    % set(gca,'XTickLength', [0 0]);
    set(gca,'XTickLabelRotation',45,'XTickLabel',{'','','Light rewarded','','Light unrewarded','','Sound rewarded','','Sound unrewarded',''})
    ylabel('Proportion approached'); %xlabel('Cue type');
    box off;
    h = gca;
    h.XRuler.TickLength = 0;
    
    % figure
    subplot(2,2,[2 4])
    gscatter([1 1 1 1],BEHAV_summary.Length.MEAN.rew_trials_light,group,'r','xo+*')
    hold on;
    gscatter([2 2 2 2],BEHAV_summary.Length.MEAN.unrew_trials_light,group,'g','xo+*')
    gscatter([3 3 3 3],BEHAV_summary.Length.MEAN.rew_trials_sound,group,'b','xo+*')
    gscatter([4 4 4 4],BEHAV_summary.Length.MEAN.unrew_trials_sound,group,'c','xo+*')
    plot([.85 1.15],[mean_rew_light2 mean_rew_light2],'k')
    plot([1.85 2.15],[mean_unrew_light2 mean_unrew_light2],'k')
    plot([2.85 3.15],[mean_rew_sound2 mean_rew_sound2],'k')
    plot([3.85 4.15],[mean_unrew_sound2 mean_unrew_sound2],'k')
    
    xlim([0 5]); ylim([0 3]); title('Trial Length');
    % set(gca,'TickLength', [0 0]); box off;
    set(gca,'XTickLabelRotation',45,'XTickLabel',{'','','Light rewarded','','Light unrewarded','','Sound rewarded','','Sound unrewarded',''})
    ylabel('Trial length (s)'); % xlabel('Cue type');
    box off;
    h = gca;
    h.XRuler.TickLength = 0;
    
    % fig = gcf;
    % fig.PaperUnits = 'inches';
    % fig.PaperPosition = [0 0 6 3];
    % print('Fig 4d - trial length','-dpng','-r0')
    %% legend
    % figure; plot([1 1 1 1],group,'xk');
    % hold on;
    % plot([1 1 1 1],group,'ok');
    % plot([1 1 1 1],group,'+k');
    % plot([1 1 1 1],group,'*k');
    % legend('show'); legend('boxoff');
    
    %% All length
    temp_R053 = length(R053_BEHAV_summary.Length.rew_trials_light);
    temp_R056 = length(R056_BEHAV_summary.Length.rew_trials_light);
    temp_R057 = length(R057_BEHAV_summary.Length.rew_trials_light);
    temp_R060 = length(R060_BEHAV_summary.Length.rew_trials_light);
    
    BEHAV_summary.Length.rew_trials_light(1:temp_R053,1) = 1;
    BEHAV_summary.Length.rew_trials_light(1:temp_R053,2) = R053_BEHAV_summary.Length.rew_trials_light;
    BEHAV_summary.Length.rew_trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.Length.rew_trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.rew_trials_light;
    BEHAV_summary.Length.rew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.Length.rew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.rew_trials_light;
    BEHAV_summary.Length.rew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.Length.rew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.rew_trials_light;
    
    BEHAV_summary.Length.unrew_trials_light(1:temp_R053,1) = 1;
    BEHAV_summary.Length.unrew_trials_light(1:temp_R053,2) = R053_BEHAV_summary.Length.unrew_trials_light;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.unrew_trials_light;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.unrew_trials_light;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.Length.unrew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.unrew_trials_light;
    
    BEHAV_summary.Length.rew_trials_sound(1:temp_R053,1) = 1;
    BEHAV_summary.Length.rew_trials_sound(1:temp_R053,2) = R053_BEHAV_summary.Length.rew_trials_sound;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.rew_trials_sound;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.rew_trials_sound;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.Length.rew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.rew_trials_sound;
    
    BEHAV_summary.Length.unrew_trials_sound(1:temp_R053,1) = 1;
    BEHAV_summary.Length.unrew_trials_sound(1:temp_R053,2) = R053_BEHAV_summary.Length.unrew_trials_sound;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.unrew_trials_sound;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.unrew_trials_sound;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.Length.unrew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.unrew_trials_sound;
    
    temp_Length_ALL = cat(1,BEHAV_summary.Length.rew_trials_light,BEHAV_summary.Length.unrew_trials_light,BEHAV_summary.Length.rew_trials_sound,BEHAV_summary.Length.unrew_trials_sound);
    temp_Length = length(BEHAV_summary.Length.rew_trials_light);
    temp_Length_rats(1:temp_Length,1) = 1;
    temp_Length_rats(temp_Length+1:temp_Length*2,1) = 2;
    temp_Length_rats(temp_Length*2+1:temp_Length*3,1) = 3;
    temp_Length_rats(temp_Length*3+1:temp_Length*4,1) = 4;
    
    BEHAV_summary.Length.ALL(:,1) = temp_Length_ALL(:,1);
    BEHAV_summary.Length.ALL(:,2) = temp_Length_rats;
    BEHAV_summary.Length.ALL(:,3) = temp_Length_ALL(:,2);
    
    %% All APP
    temp_R053 = length(R053_BEHAV_summary.APP.rew_trials_light);
    temp_R056 = length(R056_BEHAV_summary.APP.rew_trials_light);
    temp_R057 = length(R057_BEHAV_summary.APP.rew_trials_light);
    temp_R060 = length(R060_BEHAV_summary.APP.rew_trials_light);
    
    BEHAV_summary.APP.rew_trials_light(1:temp_R053,1) = 1;
    BEHAV_summary.APP.rew_trials_light(1:temp_R053,2) = R053_BEHAV_summary.APP.rew_trials_light;
    BEHAV_summary.APP.rew_trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.APP.rew_trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.rew_trials_light;
    BEHAV_summary.APP.rew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.APP.rew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.rew_trials_light;
    BEHAV_summary.APP.rew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.APP.rew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.rew_trials_light;
    
    BEHAV_summary.APP.unrew_trials_light(1:temp_R053,1) = 1;
    BEHAV_summary.APP.unrew_trials_light(1:temp_R053,2) = R053_BEHAV_summary.APP.unrew_trials_light;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.unrew_trials_light;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.unrew_trials_light;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.APP.unrew_trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.unrew_trials_light;
    
    BEHAV_summary.APP.rew_trials_sound(1:temp_R053,1) = 1;
    BEHAV_summary.APP.rew_trials_sound(1:temp_R053,2) = R053_BEHAV_summary.APP.rew_trials_sound;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.rew_trials_sound;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.rew_trials_sound;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.APP.rew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.rew_trials_sound;
    
    BEHAV_summary.APP.unrew_trials_sound(1:temp_R053,1) = 1;
    BEHAV_summary.APP.unrew_trials_sound(1:temp_R053,2) = R053_BEHAV_summary.APP.unrew_trials_sound;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.unrew_trials_sound;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.unrew_trials_sound;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
    BEHAV_summary.APP.unrew_trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.unrew_trials_sound;
    
    temp_APP_ALL = cat(1,BEHAV_summary.APP.rew_trials_light,BEHAV_summary.APP.unrew_trials_light,BEHAV_summary.APP.rew_trials_sound,BEHAV_summary.APP.unrew_trials_sound);
    temp_APP = length(BEHAV_summary.APP.rew_trials_light);
    temp_APP_rats(1:temp_APP,1) = 1;
    temp_APP_rats(temp_APP+1:temp_APP*2,1) = 2;
    temp_APP_rats(temp_APP*2+1:temp_APP*3,1) = 3;
    temp_APP_rats(temp_APP*3+1:temp_APP*4,1) = 4;
    
    for iRat = 1:temp_APP
        temp_APP_cuetype{iRat,1} = 'A';
    end
    for iRat = temp_APP+1:temp_APP*2
        temp_APP_cuetype{iRat,1} = 'B';
    end
    for iRat = temp_APP*2+1:temp_APP*3
        temp_APP_cuetype{iRat,1} = 'C';
    end
    for iRat = temp_APP*3+1:temp_APP*4
        temp_APP_cuetype{iRat,1} = 'D';
    end
    
    BEHAV_summary.APP.ALL(:,1) = temp_APP_ALL(:,1);
    BEHAV_summary.APP.ALL(:,2) = temp_APP_rats;
    BEHAV_summary.APP.ALL(:,3) = temp_APP_ALL(:,2);
    
    %% linear mixed effects model
    Length_tbl = table(BEHAV_summary.Length.ALL(:,1),BEHAV_summary.Length.ALL(:,2),BEHAV_summary.Length.ALL(:,3),'VariableNames',{'RatID','CueType','TrialLength'});
    Length_lme = fitlme(Length_tbl,'TrialLength~CueType+(1|RatID)');% +(CueType-1|RatID)'); <- matlab includes this, not sure why %lme for Trial length with cue type as fixed effect and a random intercept for each rat.
    Length_lme_reduced = fitlme(Length_tbl,'TrialLength~1+(1|RatID)');
    
    Length_comparison = compare(Length_lme_reduced,Length_lme);
    
    APP_tbl = table(BEHAV_summary.APP.ALL(:,1),temp_APP_rats,BEHAV_summary.APP.ALL(:,3),'VariableNames',{'RatID','CueType','AppProp'});
    APP_tbl.CueType = nominal(APP_tbl.CueType);
    APP_lme = fitlme(APP_tbl,'AppProp~CueType+(1|RatID)'); %+(CueType-1|RatID)');
    APP_lme_reduced = fitlme(APP_tbl,'AppProp~1+(1|RatID)');
    
    APP_comparison = compare(APP_lme_reduced,APP_lme);
    
    %% separate light and sound (ctrl) , rew and unrew exp)
    for iCue = 1:length(BEHAV_summary.APP.ALL)
        switch BEHAV_summary.APP.ALL(iCue,2)
            case {1,2}
                BEHAV_summary.APP.ALL(iCue,4) = 1;
            case {3,4}
                BEHAV_summary.APP.ALL(iCue,4) = 2;
        end
        switch BEHAV_summary.APP.ALL(iCue,2)
            case {1,3}
                BEHAV_summary.APP.ALL(iCue,5) = 1;
            case {2,4}
                BEHAV_summary.APP.ALL(iCue,5) = 2;
        end
    end
    %%
    
    APP.tbl = table(BEHAV_summary.APP.ALL(:,1),BEHAV_summary.APP.ALL(:,4),BEHAV_summary.APP.ALL(:,5),BEHAV_summary.APP.ALL(:,3),'VariableNames',{'RatID','CueIdentity','CueOutcome','AppProp'});
    APP.tbl.CueIdentity = nominal(APP.tbl.CueIdentity);
    APP.tbl.CueOutcome = nominal(APP.tbl.CueOutcome);
    APP.lme = fitlme(APP.tbl,'AppProp~CueIdentity+CueOutcome+(1|RatID)');
    APP.lme_reduced = fitlme(APP.tbl,'AppProp~CueIdentity+(1|RatID)');
    
    APP.comparison = compare(APP.lme_reduced,APP.lme);
    
    %% separate light and sound (ctrl) , rew and unrew exp)
    for iCue = 1:length(BEHAV_summary.Length.ALL)
        switch BEHAV_summary.Length.ALL(iCue,2)
            case {1,2}
                BEHAV_summary.Length.ALL(iCue,4) = 1;
            case {3,4}
                BEHAV_summary.Length.ALL(iCue,4) = 2;
        end
        switch BEHAV_summary.Length.ALL(iCue,2)
            case {1,3}
                BEHAV_summary.Length.ALL(iCue,5) = 1;
            case {2,4}
                BEHAV_summary.Length.ALL(iCue,5) = 2;
        end
    end
    %%
    
    Length.tbl = table(BEHAV_summary.Length.ALL(:,1),BEHAV_summary.Length.ALL(:,4),BEHAV_summary.Length.ALL(:,5),BEHAV_summary.Length.ALL(:,3),'VariableNames',{'RatID','CueIdentity','CueOutcome','TrialLength'});
    Length.tbl.CueIdentity = nominal(Length.tbl.CueIdentity);
    Length.tbl.CueOutcome = nominal(Length.tbl.CueOutcome);
    Length.lme = fitlme(Length.tbl,'TrialLength~CueIdentity+CueOutcome+(1|RatID)');
    Length.lme_reduced = fitlme(Length.tbl,'TrialLength~CueIdentity+(1|RatID)');
    
    Length.comparison = compare(Length.lme_reduced,Length.lme);