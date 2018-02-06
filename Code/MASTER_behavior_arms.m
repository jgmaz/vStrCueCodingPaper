% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R057\Mat2'));
Rat_id = sesh.session_id(1:4);
if strcmp(Rat_id,'R060') == 1
    block = {'TrialInfo{1,1}','TrialInfo{1,2}'};
else
    block = {'TrialInfo_block1','TrialInfo_block2'};
end
Outcome = {'unrew','rew'};    
%% rewarded vs unrewarded trials for both arms
for iArm = 1:4
    for iBlock = 1:2
        for iRew = 1:2
            count.(Outcome{iRew}).trial_num.arm{iArm} = 1;
            app.(Outcome{iRew}).trials.arm{iArm} = [];
            length.(Outcome{iRew}).trials.arm{iArm} = [];
            avg.(Outcome{iRew}).trials.arm{iArm} = [];
            SEM.(Outcome{iRew}).trials.arm{iArm} = [];
            
            for iTrial = 1:numel(metadata.(block{iBlock}).trialT)
                switch metadata.(block{iBlock}).rewarded(iTrial) %add to rewarded or unrewarded count depending on if trial was rewarded or not
                    case iRew - 1
                        if metadata.(block{iBlock}).photosensorID(iTrial) == iArm
                            app.(Outcome{iRew}).trials.arm{iArm}(count.(Outcome{iRew}).trial_num.arm{iArm}) = metadata.(block{iBlock}).summary(iTrial,3);
                            if metadata.(block{iBlock}).summary(iTrial,15) < 5
                                length.(Outcome{iRew}).trials.arm{iArm}(count.(Outcome{iRew}).trial_num.arm{iArm}) = metadata.(block{iBlock}).summary(iTrial,15);
                            else
                                length.(Outcome{iRew}).trials.arm{iArm}(count.(Outcome{iRew}).trial_num.arm{iArm}) = NaN;
                            end
                            count.(Outcome{iRew}).trial_num.arm{iArm} = count.(Outcome{iRew}).trial_num.arm{iArm} + 1;
                        end
                end
            end
            
        end
    end
end
        
%% summary

count.(Outcome{iRew}).trials.arm{iArm} = length(app.(Outcome{iRew}).trials.arm{iArm});
prop_app.(Outcome{iRew}).trials.arm{iArm} = sum(app.(Outcome{iRew}).trials.arm{iArm})/count.(Outcome{iRew}).trials.arm{iArm};
avg_length.(Outcome{iRew}).trials.arm{iArm} = nanmean(length.(Outcome{iRew}).trials.arm{iArm});
SEM_length.(Outcome{iRew}).trials.arm{iArm} = nanstd(length.(Outcome{iRew}).trials.arm{iArm}/sqrt(numel(length.(Outcome{iRew}).trials.arm{iArm})-sum(isnan(length.(Outcome{iRew}).trials.arm{iArm}))));

count_(Outcome{iRew}).trials.arm{iArm} = length(app_(Outcome{iRew}).trials.arm{iArm});
prop_app_(Outcome{iRew}).trials.arm{iArm} = sum(app_(Outcome{iRew}).trials.arm{iArm})/count_(Outcome{iRew}).trials.arm{iArm};
avg_length_(Outcome{iRew}).trials.arm{iArm} = nanmean(length_(Outcome{iRew}).trials.arm{iArm});
SEM_length_(Outcome{iRew}).trials.arm{iArm} = nanstd(length_(Outcome{iRew}).trials.arm{iArm}/sqrt(numel(length_(Outcome{iRew}).trials.arm{iArm})-sum(isnan(length_(Outcome{iRew}).trials.arm{iArm}))));

count.(Outcome{iRew}).trials.arm{iArm} = length(app.(Outcome{iRew}).trials.arm{iArm});
prop_app.(Outcome{iRew}).trials.arm{iArm} = sum(app.(Outcome{iRew}).trials.arm{iArm})/count.(Outcome{iRew}).trials.arm{iArm};
avg_length.(Outcome{iRew}).trials.arm{iArm} = nanmean(length.(Outcome{iRew}).trials.arm{iArm});
SEM_length.(Outcome{iRew}).trials.arm{iArm} = nanstd(length.(Outcome{iRew}).trials.arm{iArm}/sqrt(numel(length.(Outcome{iRew}).trials.arm{iArm})-sum(isnan(length.(Outcome{iRew}).trials.arm{iArm}))));

count_(Outcome{iRew}).trials.arm{iArm} = length(app_(Outcome{iRew}).trials.arm{iArm});
prop_app_(Outcome{iRew}).trials.arm{iArm} = sum(app_(Outcome{iRew}).trials.arm{iArm})/count_(Outcome{iRew}).trials.arm{iArm};
avg_length_(Outcome{iRew}).trials.arm{iArm} = nanmean(length_(Outcome{iRew}).trials.arm{iArm});
SEM_length_(Outcome{iRew}).trials.arm{iArm} = nanstd(length_(Outcome{iRew}).trials.arm{iArm}/sqrt(numel(length_(Outcome{iRew}).trials.arm{iArm})-sum(isnan(length_(Outcome{iRew}).trials.arm{iArm}))));

%% output
switch sesh.arm_order
    case 1
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_light = avg_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_light = avg_length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound = avg_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound = avg_length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Length.SEM.(Outcome{iRew}).trials_light = SEM_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_light = SEM_length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_sound = SEM_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_sound = SEM_length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Length.ALL.(Outcome{iRew}).trials_light = length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_light = length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_sound = length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_sound = length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_light = count.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_light = count_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_sound = count.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_sound = count_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_light = prop_app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_light = prop_app_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound = prop_app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound = prop_app_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_light = app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_light = app_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_sound = app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_sound = app_(Outcome{iRew}).trials.arm{iArm};
        
    case 2
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_light = avg_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_light = avg_length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound = avg_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound = avg_length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Length.SEM.(Outcome{iRew}).trials_light = SEM_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_light = SEM_length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_sound = SEM_length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.SEM.(Outcome{iRew}).trials_sound = SEM_length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Length.ALL.(Outcome{iRew}).trials_light = length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_light = length_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_sound = length.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Length.ALL.(Outcome{iRew}).trials_sound = length_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_light = count.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_light = count_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_sound = count.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.COUNT.(Outcome{iRew}).trials_sound = count_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_light = prop_app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_light = prop_app_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound = prop_app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound = prop_app_(Outcome{iRew}).trials.arm{iArm};
        
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_light = app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_light = app_(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_sound = app.(Outcome{iRew}).trials.arm{iArm};
        BEHAV.Approach.ALL.(Outcome{iRew}).trials_sound = app_(Outcome{iRew}).trials.arm{iArm};
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
end
clear
%% group
% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R053\Mat\BEHAV\BEHAV'));

counter = 1;
mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(num2str(kk));
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
    BEHAV_summary.APP.(Outcome{iRew}).trials_light(kk) =  BEHAV.Approach.PROP.(Outcome{iRew}).trials_light;
    BEHAV_summary.APP.(Outcome{iRew}).trials_light(kk) =  BEHAV.Approach.PROP.(Outcome{iRew}).trials_light;
    BEHAV_summary.APP.(Outcome{iRew}).trials_sound(kk) =  BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound;
    BEHAV_summary.APP.(Outcome{iRew}).trials_sound(kk) =  BEHAV.Approach.PROP.(Outcome{iRew}).trials_sound;
    BEHAV_summary.Length.(Outcome{iRew}).trials_light(kk) =  BEHAV.Length.MEAN.(Outcome{iRew}).trials_light;
    BEHAV_summary.Length.(Outcome{iRew}).trials_light(kk) =  BEHAV.Length.MEAN.(Outcome{iRew}).trials_light;
    BEHAV_summary.Length.(Outcome{iRew}).trials_sound(kk) =  BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound;
    BEHAV_summary.Length.(Outcome{iRew}).trials_sound(kk) =  BEHAV.Length.MEAN.(Outcome{iRew}).trials_sound;
end

%% group summary

BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light = mean(BEHAV_summary.APP.(Outcome{iRew}).trials_light);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light = mean(BEHAV_summary.APP.(Outcome{iRew}).trials_light);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound = mean(BEHAV_summary.APP.(Outcome{iRew}).trials_sound);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound = mean(BEHAV_summary.APP.(Outcome{iRew}).trials_sound);

BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light = mean(BEHAV_summary.Length.(Outcome{iRew}).trials_light);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light = mean(BEHAV_summary.Length.(Outcome{iRew}).trials_light);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound = mean(BEHAV_summary.Length.(Outcome{iRew}).trials_sound);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound = mean(BEHAV_summary.Length.(Outcome{iRew}).trials_sound);

BEHAV_summary.APP.SEM.(Outcome{iRew}).trials_light = std(BEHAV_summary.APP.(Outcome{iRew}).trials_light)/sqrt(numel(BEHAV_summary.APP.(Outcome{iRew}).trials_light));
BEHAV_summary.APP.SEM.(Outcome{iRew}).trials_light = std(BEHAV_summary.APP.(Outcome{iRew}).trials_light)/sqrt(numel(BEHAV_summary.APP.(Outcome{iRew}).trials_light));
BEHAV_summary.APP.SEM.(Outcome{iRew}).trials_sound = std(BEHAV_summary.APP.(Outcome{iRew}).trials_sound)/sqrt(numel(BEHAV_summary.APP.(Outcome{iRew}).trials_sound));
BEHAV_summary.APP.SEM.(Outcome{iRew}).trials_sound = std(BEHAV_summary.APP.(Outcome{iRew}).trials_sound)/sqrt(numel(BEHAV_summary.APP.(Outcome{iRew}).trials_sound));

BEHAV_summary.Length.SEM.(Outcome{iRew}).trials_light = std(BEHAV_summary.Length.(Outcome{iRew}).trials_light)/sqrt(numel(BEHAV_summary.Length.(Outcome{iRew}).trials_light));
BEHAV_summary.Length.SEM.(Outcome{iRew}).trials_light = std(BEHAV_summary.Length.(Outcome{iRew}).trials_light)/sqrt(numel(BEHAV_summary.Length.(Outcome{iRew}).trials_light));
BEHAV_summary.Length.SEM.(Outcome{iRew}).trials_sound = std(BEHAV_summary.Length.(Outcome{iRew}).trials_sound)/sqrt(numel(BEHAV_summary.Length.(Outcome{iRew}).trials_sound));
BEHAV_summary.Length.SEM.(Outcome{iRew}).trials_sound = std(BEHAV_summary.Length.(Outcome{iRew}).trials_sound)/sqrt(numel(BEHAV_summary.Length.(Outcome{iRew}).trials_sound));

%%
BEHAV_summary.APP.ALL(:,1) = BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.ALL(:,2) = BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.ALL(:,3) = BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.ALL(:,4) = BEHAV_summary.APP.(Outcome{iRew}).trials_sound;

BEHAV_summary.Length.ALL(:,1) = BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.ALL(:,2) = BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.ALL(:,3) = BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.ALL(:,4) = BEHAV_summary.Length.(Outcome{iRew}).trials_sound;

%%

[p_APP,tbl_APP,stats_APP] = anova1(BEHAV_summary.APP.ALL);
[p_Length,tbl_Length,stats_Length] = anova1(BEHAV_summary.Length.ALL);


[h_rvr p_rvr] = ttest(BEHAV_summary.APP.(Outcome{iRew}).trials_sound,BEHAV_summary.APP.(Outcome{iRew}).trials_light);
[h_uvu p_uvu] = ttest(BEHAV_summary.APP.(Outcome{iRew}).trials_sound,BEHAV_summary.APP.(Outcome{iRew}).trials_light);
[h_rvu1 p_rvu1] = ttest(BEHAV_summary.APP.(Outcome{iRew}).trials_light,BEHAV_summary.APP.(Outcome{iRew}).trials_light);
[h_rvu2 p_rvu2] = ttest(BEHAV_summary.APP.(Outcome{iRew}).trials_sound,BEHAV_summary.APP.(Outcome{iRew}).trials_sound);

[lh_rvr lp_rvr] = ttest(BEHAV_summary.Length.(Outcome{iRew}).trials_sound,BEHAV_summary.Length.(Outcome{iRew}).trials_light);
[lh_uvu lp_uvu] = ttest(BEHAV_summary.Length.(Outcome{iRew}).trials_sound,BEHAV_summary.Length.(Outcome{iRew}).trials_light);
[lh_rvu1 lp_rvu1] = ttest(BEHAV_summary.Length.(Outcome{iRew}).trials_light,BEHAV_summary.Length.(Outcome{iRew}).trials_light);
[lh_rvu2 lp_rvu2] = ttest(BEHAV_summary.Length.(Outcome{iRew}).trials_sound,BEHAV_summary.Length.(Outcome{iRew}).trials_sound);

%% all rats together
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light = cat(2,R053_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light, ...
    R056_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,R057_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,R060_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light = cat(2,R053_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light, ...
    R056_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,R057_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,R060_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound = cat(2,R053_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound, ...
    R056_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,R057_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,R060_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound);
BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound = cat(2,R053_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound, ...
    R056_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,R057_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,R060_BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound);

BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light = cat(2,R053_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light, ...
    R056_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,R057_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,R060_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light = cat(2,R053_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light, ...
    R056_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,R057_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,R060_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound = cat(2,R053_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound, ...
    R056_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,R057_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,R060_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound);
BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound = cat(2,R053_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound, ...
    R056_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,R057_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,R060_BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound);

%%
mean.(Outcome{iRew}).light = mean(BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light);
mean_(Outcome{iRew}).light = mean(BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light);
mean.(Outcome{iRew}).sound = mean(BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound);
mean_(Outcome{iRew}).sound = mean(BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound);

mean.(Outcome{iRew}).light2 = mean(BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light);
mean_(Outcome{iRew}).light2 = mean(BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light);
mean.(Outcome{iRew}).sound2 = mean(BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound);
mean_(Outcome{iRew}).sound2 = mean(BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound);
%%
group = [1 2 3 4];
figure
subplot(2,2,[1 3])
gscatter([1 1 1 1],BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,group,'r','xo+*')
hold on;
gscatter([2 2 2 2],BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_light,group,'g','xo+*')
gscatter([3 3 3 3],BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,group,'b','xo+*')
gscatter([4 4 4 4],BEHAV_summary.APP.MEAN.(Outcome{iRew}).trials_sound,group,'c','xo+*')
plot([.85 1.15],[mean.(Outcome{iRew}).light mean.(Outcome{iRew}).light],'k')
plot([1.85 2.15],[mean_(Outcome{iRew}).light mean_(Outcome{iRew}).light],'k')
plot([2.85 3.15],[mean.(Outcome{iRew}).sound mean.(Outcome{iRew}).sound],'k')
plot([3.85 4.15],[mean_(Outcome{iRew}).sound mean_(Outcome{iRew}).sound],'k')

xlim([0 5]); ylim([0 1]); title('Proportion Approached');
% set(gca,'XTickLength', [0 0]); 
set(gca,'XTickLabelRotation',45,'XTickLabel',{'','','Light rewarded','','Light unrewarded','','Sound rewarded','','Sound unrewarded',''})
ylabel('Proportion approached'); %xlabel('Cue type');
box off;
h = gca;
h.XRuler.TickLength = 0;  
    
% figure
subplot(2,2,[2 4])
gscatter([1 1 1 1],BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,group,'r','xo+*')
hold on;
gscatter([2 2 2 2],BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_light,group,'g','xo+*')
gscatter([3 3 3 3],BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,group,'b','xo+*')
gscatter([4 4 4 4],BEHAV_summary.Length.MEAN.(Outcome{iRew}).trials_sound,group,'c','xo+*')
plot([.85 1.15],[mean.(Outcome{iRew}).light2 mean.(Outcome{iRew}).light2],'k')
plot([1.85 2.15],[mean_(Outcome{iRew}).light2 mean_(Outcome{iRew}).light2],'k')
plot([2.85 3.15],[mean.(Outcome{iRew}).sound2 mean.(Outcome{iRew}).sound2],'k')
plot([3.85 4.15],[mean_(Outcome{iRew}).sound2 mean_(Outcome{iRew}).sound2],'k')

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
temp_R053 = length(R053_BEHAV_summary.Length.(Outcome{iRew}).trials_light);
temp_R056 = length(R056_BEHAV_summary.Length.(Outcome{iRew}).trials_light);
temp_R057 = length(R057_BEHAV_summary.Length.(Outcome{iRew}).trials_light);
temp_R060 = length(R060_BEHAV_summary.Length.(Outcome{iRew}).trials_light);

BEHAV_summary.Length.(Outcome{iRew}).trials_light(1:temp_R053,1) = 1;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(1:temp_R053,2) = R053_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.(Outcome{iRew}).trials_light;

BEHAV_summary.Length.(Outcome{iRew}).trials_light(1:temp_R053,1) = 1;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(1:temp_R053,2) = R053_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.(Outcome{iRew}).trials_light;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.Length.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.(Outcome{iRew}).trials_light;

BEHAV_summary.Length.(Outcome{iRew}).trials_sound(1:temp_R053,1) = 1;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(1:temp_R053,2) = R053_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;

BEHAV_summary.Length.(Outcome{iRew}).trials_sound(1:temp_R053,1) = 1;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(1:temp_R053,2) = R053_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.Length.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.Length.(Outcome{iRew}).trials_sound;

temp_Length_ALL = cat(1,BEHAV_summary.Length.(Outcome{iRew}).trials_light,BEHAV_summary.Length.(Outcome{iRew}).trials_light,BEHAV_summary.Length.(Outcome{iRew}).trials_sound,BEHAV_summary.Length.(Outcome{iRew}).trials_sound);
temp_Length = length(BEHAV_summary.Length.(Outcome{iRew}).trials_light);
temp_Length_rats(1:temp_Length,1) = 1;
temp_Length_rats(temp_Length+1:temp_Length*2,1) = 2;
temp_Length_rats(temp_Length*2+1:temp_Length*3,1) = 3;
temp_Length_rats(temp_Length*3+1:temp_Length*4,1) = 4;

BEHAV_summary.Length.ALL(:,1) = temp_Length_ALL(:,1);
BEHAV_summary.Length.ALL(:,2) = temp_Length_rats;
BEHAV_summary.Length.ALL(:,3) = temp_Length_ALL(:,2);

%% All APP
temp_R053 = length(R053_BEHAV_summary.APP.(Outcome{iRew}).trials_light);
temp_R056 = length(R056_BEHAV_summary.APP.(Outcome{iRew}).trials_light);
temp_R057 = length(R057_BEHAV_summary.APP.(Outcome{iRew}).trials_light);
temp_R060 = length(R060_BEHAV_summary.APP.(Outcome{iRew}).trials_light);

BEHAV_summary.APP.(Outcome{iRew}).trials_light(1:temp_R053,1) = 1;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(1:temp_R053,2) = R053_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.(Outcome{iRew}).trials_light;

BEHAV_summary.APP.(Outcome{iRew}).trials_light(1:temp_R053,1) = 1;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(1:temp_R053,2) = R053_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.(Outcome{iRew}).trials_light;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.APP.(Outcome{iRew}).trials_light(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.(Outcome{iRew}).trials_light;

BEHAV_summary.APP.(Outcome{iRew}).trials_sound(1:temp_R053,1) = 1;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(1:temp_R053,2) = R053_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;

BEHAV_summary.APP.(Outcome{iRew}).trials_sound(1:temp_R053,1) = 1;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(1:temp_R053,2) = R053_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,1) = 2;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+1:temp_R053+temp_R056,2) = R056_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,1) = 3;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+1:temp_R053+temp_R056+temp_R057,2) = R057_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,1) = 4;
BEHAV_summary.APP.(Outcome{iRew}).trials_sound(temp_R053+temp_R056+temp_R057+1:temp_R053+temp_R056+temp_R057+temp_R060,2) = R060_BEHAV_summary.APP.(Outcome{iRew}).trials_sound;

temp_APP_ALL = cat(1,BEHAV_summary.APP.(Outcome{iRew}).trials_light,BEHAV_summary.APP.(Outcome{iRew}).trials_light,BEHAV_summary.APP.(Outcome{iRew}).trials_sound,BEHAV_summary.APP.(Outcome{iRew}).trials_sound);
temp_APP = length(BEHAV_summary.APP.(Outcome{iRew}).trials_light);
temp_APP_rats(1:temp_APP,1) = 1;
temp_APP_rats(temp_APP+1:temp_APP*2,1) = 2;
temp_APP_rats(temp_APP*2+1:temp_APP*3,1) = 3;
temp_APP_rats(temp_APP*3+1:temp_APP*4,1) = 4;

for iRat = 1:temp_APP
    temp_APP_rats{iRat,1} = 'A';
end
for iRat = temp_APP+1:temp_APP*2
    temp_APP_rats{iRat,1} = 'B';
end
for iRat = temp_APP*2+1:temp_APP*3
    temp_APP_rats{iRat,1} = 'C';
end
for iRat = temp_APP*3+1:temp_APP*4
    temp_APP_rats{iRat,1} = 'D';
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