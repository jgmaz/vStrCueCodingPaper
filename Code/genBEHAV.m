function BEHAV = genBEHAV(directory,destination)
% function BEHAV = genBEHAV(directory,destination)
%
%
% INPUTS:
%
% OUTPUTS:
%
f = filesep;

iSesh = 1;
for ii = 1:4 % work through each rat in the dataset
    switch ii
        case 1
            rat_id = 'R053';
            day = {'-2014-11-04'; '-2014-11-05'; '-2014-11-06'; '-2014-11-07'; '-2014-11-08'; ...
                '-2014-11-10'; '-2014-11-11'; '-2014-11-12'; '-2014-11-13'; '-2014-11-14'; ...
                '-2014-11-15'; '-2014-11-16'};
            block_order_list = [1 1 2 1 2 2 1 2 2 1 1 2];
        case 2
            rat_id = 'R056';
            day = {'-2015-05-16'; '-2015-05-18'; '-2015-05-19'; '-2015-05-21'; '-2015-05-22'; ...
                '-2015-05-23'; '-2015-05-29'; '-2015-05-31'; '-2015-06-01'; '-2015-06-02'; ...
                '-2015-06-03'; '-2015-06-04'; '-2015-06-05'; '-2015-06-08'; '-2015-06-09'};
            block_order_list = [1 1 2 2 2 1 1 2 2 1 2 1 2 1 2];
        case 3
            rat_id = 'R057';
            day = {'-2015-02-14'; '-2015-02-15'; '-2015-02-16'; '-2015-02-17'; '-2015-02-18'; ...
                '-2015-02-21'; '-2015-02-22'; '-2015-02-23'; '-2015-02-24'; '-2015-02-25'; ...
                '-2015-02-26'; '-2015-02-27'};
            block_order_list = [2 1 2 1 2 1 2 1 2 1 2 1];
        case 4
            rat_id = 'R060';
            day = {'-2014-12-12'; '-2014-12-13'; '-2014-12-14'; '-2014-12-15'; ...
                '-2014-12-17'; '-2014-12-20'; '-2014-12-23'; '-2014-12-24'; ...
                '-2014-12-26'; '-2014-12-27'; '-2014-12-28'; '-2014-12-29'; '-2014-12-30'; '-2014-12-31'; ...
                '-2015-01-01'; '-2015-01-02'; '-2015-01-03'; '-2015-01-04'};
            block_order_list = [1 2 1 2 1 2 1 2 1 1 2 1 1 2 1 1 2 1 2 1 2];
    end
    for kj = 1:length(day)
        
        %% input session info
        sesh.session_id = cat(2,rat_id,day{kj}); % current day for analysis
        
        sesh.block_order = block_order_list(kj); % 1 if light block came first, 2 if sound block came first
        disp(cat(2,'loading session ',sesh.session_id));
        cd(cat(2,directory,sesh.session_id(1:4),f,sesh.session_id))
        load(strcat(sesh.session_id,'_metadata.mat'))
        
        rew_trial_num_block1 = 1;
        unrew_trial_num_block1 = 1;
        rew_trial_num_block2 = 1;
        unrew_trial_num_block2 = 1;
        
        app_rew_trials_block1 = [];
        app_unrew_trials_block1 = [];
        app_rew_trials_block2 = [];
        app_unrew_trials_block2 = [];
        
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
                    rew_trial_num_block1 = rew_trial_num_block1 + 1;
                case 0
                    app_unrew_trials_block1(unrew_trial_num_block1) = metadata.TrialInfo{1,1}.summary(ik,3);
                    unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
            end
        end
        
        for ip = 1:length(metadata.TrialInfo{1,2}.trialT)
            switch metadata.TrialInfo{1,2}.rewarded(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
                case 1
                    app_rew_trials_block2(rew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,3);
                    rew_trial_num_block2 = rew_trial_num_block2 + 1;
                case 0
                    app_unrew_trials_block2(unrew_trial_num_block2) = metadata.TrialInfo{1,2}.summary(ip,3);
                    unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
            end
        end
        %% summary
        
        count_rew_trials_block1 = length(app_rew_trials_block1);
        prop_app_rew_trials_block1 = sum(app_rew_trials_block1)/count_rew_trials_block1;
        
        count_unrew_trials_block1 = length(app_unrew_trials_block1);
        prop_app_unrew_trials_block1 = sum(app_unrew_trials_block1)/count_unrew_trials_block1;
        
        count_rew_trials_block2 = length(app_rew_trials_block2);
        prop_app_rew_trials_block2 = sum(app_rew_trials_block2)/count_rew_trials_block2;
        
        count_unrew_trials_block2 = length(app_unrew_trials_block2);
        prop_app_unrew_trials_block2 = sum(app_unrew_trials_block2)/count_unrew_trials_block2;
        
        %% output
        switch sesh.block_order
            case 1
                BEHAV.Approach.PROP.rew_trials_light = prop_app_rew_trials_block1;
                BEHAV.Approach.PROP.unrew_trials_light = prop_app_unrew_trials_block1;
                BEHAV.Approach.PROP.rew_trials_sound = prop_app_rew_trials_block2;
                BEHAV.Approach.PROP.unrew_trials_sound = prop_app_unrew_trials_block2;
                
            case 2
                BEHAV.Approach.PROP.rew_trials_light = prop_app_rew_trials_block2;
                BEHAV.Approach.PROP.unrew_trials_light = prop_app_unrew_trials_block2;
                BEHAV.Approach.PROP.rew_trials_sound = prop_app_rew_trials_block1;
                BEHAV.Approach.PROP.unrew_trials_sound = prop_app_unrew_trials_block1;
                
        end
        
        BEHAV_summary.APP.rew_trials_light(iSesh,1) =  BEHAV.Approach.PROP.rew_trials_light;
        BEHAV_summary.APP.unrew_trials_light(iSesh,1) =  BEHAV.Approach.PROP.unrew_trials_light;
        BEHAV_summary.APP.rew_trials_sound(iSesh,1) =  BEHAV.Approach.PROP.rew_trials_sound;
        BEHAV_summary.APP.unrew_trials_sound(iSesh,1) =  BEHAV.Approach.PROP.unrew_trials_sound;
        BEHAV_summary.APP.ratID(iSesh,1) = ii;
        
        clearvars -except f BEHAV_summary iSesh day kj ii rat_id directory destination block_order_list
        iSesh = iSesh + 1;
    end
end

%%
BEHAV_summary.APP.ALL(:,1) = cat(1,BEHAV_summary.APP.ratID,BEHAV_summary.APP.ratID,BEHAV_summary.APP.ratID,BEHAV_summary.APP.ratID);
BEHAV_summary.APP.ALL(:,3) = cat(1,BEHAV_summary.APP.rew_trials_light,BEHAV_summary.APP.unrew_trials_light,BEHAV_summary.APP.rew_trials_sound,BEHAV_summary.APP.unrew_trials_sound);

temp_APP = length(BEHAV_summary.APP.rew_trials_light);
temp_APP_cue(1:temp_APP,1) = 1;
temp_APP_cue(temp_APP+1:temp_APP*2,1) = 2;
temp_APP_cue(temp_APP*2+1:temp_APP*3,1) = 3;
temp_APP_cue(temp_APP*3+1:temp_APP*4,1) = 4;
BEHAV_summary.APP.ALL(:,2) = temp_APP_cue;

%% separate light and sound (ctrl) , rew and unrew (exp)
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
for iRat = 1:4
    idx = find(BEHAV_summary.APP.ALL(:,1) == iRat & BEHAV_summary.APP.ALL(:,2) == 1);
    BEHAV_summary.APP.MEAN.rew_trials_light(iRat) = mean(BEHAV_summary.APP.ALL(idx,3));
    idx = find(BEHAV_summary.APP.ALL(:,1) == iRat & BEHAV_summary.APP.ALL(:,2) == 2);
    BEHAV_summary.APP.MEAN.unrew_trials_light(iRat) = mean(BEHAV_summary.APP.ALL(idx,3));
    idx = find(BEHAV_summary.APP.ALL(:,1) == iRat & BEHAV_summary.APP.ALL(:,2) == 3);
    BEHAV_summary.APP.MEAN.rew_trials_sound(iRat) = mean(BEHAV_summary.APP.ALL(idx,3));
    idx = find(BEHAV_summary.APP.ALL(:,1) == iRat & BEHAV_summary.APP.ALL(:,2) == 4);
    BEHAV_summary.APP.MEAN.unrew_trials_sound(iRat) = mean(BEHAV_summary.APP.ALL(idx,3));
end

%%
BEHAV_summary.mean_rew_light = mean(BEHAV_summary.APP.MEAN.rew_trials_light);
BEHAV_summary.mean_unrew_light = mean(BEHAV_summary.APP.MEAN.unrew_trials_light);
BEHAV_summary.mean_rew_sound = mean(BEHAV_summary.APP.MEAN.rew_trials_sound);
BEHAV_summary.mean_unrew_sound = mean(BEHAV_summary.APP.MEAN.unrew_trials_sound);

%%

APP.tbl = table(BEHAV_summary.APP.ALL(:,1),BEHAV_summary.APP.ALL(:,4),BEHAV_summary.APP.ALL(:,5),BEHAV_summary.APP.ALL(:,3),'VariableNames',{'RatID','CueIdentity','CueOutcome','AppProp'});
APP.tbl.CueIdentity = nominal(APP.tbl.CueIdentity);
APP.tbl.CueOutcome = nominal(APP.tbl.CueOutcome);
APP.lme = fitlme(APP.tbl,'AppProp~CueIdentity+CueOutcome+(1|RatID)');
APP.lme_reduced = fitlme(APP.tbl,'AppProp~CueIdentity+(1|RatID)');

APP.comparison = compare(APP.lme_reduced,APP.lme);

save(strcat(destination,'Behavior_summary.mat'),'BEHAV_summary','APP')

end