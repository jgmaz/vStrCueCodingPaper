
    counter = 1;
    mat_files = dir('*.mat');
    for kk = 1:length(dir('*.mat'))
        disp(num2str(kk));
        load(mat_files(kk).name);
        mat_overview.fname{kk} = mat_files(kk).name;
    
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
    
        new_v_old = strcmp(mat_overview.fname{kk}(1:4),'R060');   
            switch new_v_old
                case 0
                    for ik = 1:length(metadata.TrialInfo_block1.trialT)
        switch metadata.TrialInfo_block1.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                app_rew_trials_block1(rew_trial_num_block1) = metadata.TrialInfo_block1.summary(ik,3);
                rew_trial_num_block1 = rew_trial_num_block1 + 1;
            case 0
                app_unrew_trials_block1(unrew_trial_num_block1) = metadata.TrialInfo_block1.summary(ik,3);
                unrew_trial_num_block1 = unrew_trial_num_block1 + 1;
        end
    end
    
    for ip = 1:length(metadata.TrialInfo_block2.trialT)
        switch metadata.TrialInfo_block2.rewarded(ip) %add to rewarded or unrewarded count depending on if trial was rewarded or not
            case 1
                app_rew_trials_block2(rew_trial_num_block2) = metadata.TrialInfo_block2.summary(ip,3);
                rew_trial_num_block2 = rew_trial_num_block2 + 1;
            case 0
                app_unrew_trials_block2(unrew_trial_num_block2) = metadata.TrialInfo_block2.summary(ip,3);
                unrew_trial_num_block2 = unrew_trial_num_block2 + 1;
        end
    end
                    
                case 1                   
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
    
    BEHAV.Approach.COUNT.rew_trials_block1 = count_rew_trials_block1;
            BEHAV.App_OG.COUNT.unrew_trials_block1 = count_unrew_trials_block1;
            BEHAV.App_OG.COUNT.rew_trials_block2 = count_rew_trials_block2;
            BEHAV.App_OG.COUNT.unrew_trials_block2 = count_unrew_trials_block2;
            
            BEHAV.App_OG.PROP.rew_trials_block1 = prop_app_rew_trials_block1;
            BEHAV.App_OG.PROP.unrew_trials_block1 = prop_app_unrew_trials_block1;
            BEHAV.App_OG.PROP.rew_trials_block2 = prop_app_rew_trials_block2;
            BEHAV.App_OG.PROP.unrew_trials_block2 = prop_app_unrew_trials_block2;
            
            BEHAV.App_OG.ALL.rew_trials_block1 = app_rew_trials_block1;
            BEHAV.App_OG.ALL.unrew_trials_block1 = app_unrew_trials_block1;
            BEHAV.App_OG.ALL.rew_trials_block2 = app_rew_trials_block2;
            BEHAV.App_OG.ALL.unrew_trials_block2 = app_unrew_trials_block2;
    
    BEHAVS.Summary.(cat(2,'sesh',num2str(kk))) = BEHAV.Approach;
    BEHAVS.OGBlock.(cat(2,'sesh',num2str(kk))) = BEHAV.App_OG;
    
     clearvars -except kk BEHAVS mat_files
    end

%%
cue_type = {'rew_trials_light' 'unrew_trials_light' 'rew_trials_sound' 'unrew_trials_sound'};
cue_type = {'rew_trials_block1' 'unrew_trials_block1' 'rew_trials_block2' 'unrew_trials_block2'};
Type = {'Half' 'Third'};
for iSesh = 1:57%length(BEHAVS.Summary)
    for iCue = 1:length(cue_type)
        n_trials = length(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue}));
        half_trials = floor(n_trials/2);
        half_2_trials = n_trials - half_trials;
        third_trials = round(n_trials/3);
        third_2_trials = n_trials - third_trials; 
        BEHAV.Total.Half.(cue_type{iCue})(iSesh) = half_trials;
        BEHAV.Total.Third.(cue_type{iCue})(iSesh) = third_trials;
        for iHalve = 1:2
            switch iHalve
                case 1
    BEHAV.Count.Half.(cue_type{iCue})(iSesh,iHalve) = sum(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue})(1:half_trials));
    BEHAV.Prop.Half.(cue_type{iCue})(iSesh,iHalve) = BEHAV.Count.Half.(cue_type{iCue})(iSesh,iHalve)/half_trials;
                case 2
                        BEHAV.Count.Half.(cue_type{iCue})(iSesh,iHalve) = sum(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue})(half_2_trials+1:end));
    BEHAV.Prop.Half.(cue_type{iCue})(iSesh,iHalve) = BEHAV.Count.Half.(cue_type{iCue})(iSesh,iHalve)/half_trials;
            end
        end
        for iThird = 1:3
            switch iThird
                case 1
    BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird) = sum(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue})(1:third_trials));
    BEHAV.Prop.Third.(cue_type{iCue})(iSesh,iThird) = BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird)/third_trials;
                case 2
                        BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird) = sum(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue})(third_trials+1:third_2_trials));
    BEHAV.Prop.Third.(cue_type{iCue})(iSesh,iThird) = BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird)/length(third_trials+1:third_2_trials);
                    case 3
                        BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird) = sum(BEHAVS.Summary.(cat(2,'sesh',num2str(iSesh))).ALL.(cue_type{iCue})(third_2_trials+1:end));
    BEHAV.Prop.Third.(cue_type{iCue})(iSesh,iThird) = BEHAV.Count.Third.(cue_type{iCue})(iSesh,iThird)/third_trials;
            end
        end
        BEHAV.Count.Half.Diff.(cue_type{iCue})(iSesh) = BEHAV.Count.Half.(cue_type{iCue})(iSesh,2) - BEHAV.Count.Half.(cue_type{iCue})(iSesh,1);
        BEHAV.Prop.Half.Diff.(cue_type{iCue})(iSesh) = BEHAV.Prop.Half.(cue_type{iCue})(iSesh,2) - BEHAV.Prop.Half.(cue_type{iCue})(iSesh,1);
        BEHAV.Count.Third.Diff.(cue_type{iCue})(iSesh) = BEHAV.Count.Third.(cue_type{iCue})(iSesh,3) - BEHAV.Count.Third.(cue_type{iCue})(iSesh,1);
        BEHAV.Prop.Third.Diff.(cue_type{iCue})(iSesh) = BEHAV.Prop.Third.(cue_type{iCue})(iSesh,3) - BEHAV.Prop.Third.(cue_type{iCue})(iSesh,1);
    end
    
    for iType = 1:2
    BEHAV.Corr.RAW.(Type{iType})(iSesh,1) = BEHAV.Count.(Type{iType}).(cue_type{1})(iSesh,1) + (BEHAV.Total.(Type{iType}).(cue_type{2})(iSesh) - BEHAV.Count.(Type{iType}).(cue_type{2})(iSesh,1));
    BEHAV.Corr.Prop.(Type{iType})(iSesh,1) = BEHAV.Corr.RAW.(Type{iType})(iSesh,1) / (BEHAV.Total.(Type{iType}).(cue_type{1})(iSesh) + BEHAV.Total.(Type{iType}).(cue_type{2})(iSesh));
    BEHAV.Corr.RAW.(Type{iType})(iSesh,2) = BEHAV.Count.(Type{iType}).(cue_type{1})(iSesh,2) + (BEHAV.Total.(Type{iType}).(cue_type{2})(iSesh) - BEHAV.Count.(Type{iType}).(cue_type{2})(iSesh,2));
    BEHAV.Corr.Prop.(Type{iType})(iSesh,2) = BEHAV.Corr.RAW.(Type{iType})(iSesh,2) / (BEHAV.Total.(Type{iType}).(cue_type{1})(iSesh) + BEHAV.Total.(Type{iType}).(cue_type{2})(iSesh));
    BEHAV.Corr.RAW.(Type{iType})(iSesh,3) = BEHAV.Count.(Type{iType}).(cue_type{3})(iSesh,1) + (BEHAV.Total.(Type{iType}).(cue_type{4})(iSesh) - BEHAV.Count.(Type{iType}).(cue_type{4})(iSesh,1));
    BEHAV.Corr.Prop.(Type{iType})(iSesh,3) = BEHAV.Corr.RAW.(Type{iType})(iSesh,3) / (BEHAV.Total.(Type{iType}).(cue_type{3})(iSesh) + BEHAV.Total.(Type{iType}).(cue_type{4})(iSesh));
    BEHAV.Corr.RAW.(Type{iType})(iSesh,4) = BEHAV.Count.(Type{iType}).(cue_type{3})(iSesh,2) + (BEHAV.Total.(Type{iType}).(cue_type{4})(iSesh) - BEHAV.Count.(Type{iType}).(cue_type{4})(iSesh,2));
    BEHAV.Corr.Prop.(Type{iType})(iSesh,4) = BEHAV.Corr.RAW.(Type{iType})(iSesh,4) / (BEHAV.Total.(Type{iType}).(cue_type{3})(iSesh) + BEHAV.Total.(Type{iType}).(cue_type{4})(iSesh));
    end
end
        
%%
% figure
% iPlot = 1;
Type = {'Half' 'Third'};
cueT = {'Block 1 +' 'Block 1 -' 'Block 2 +' 'Block 2 -'};
last = [2 3];
for iType = 1:2
    figure
    iPlot = 1;
    for iStyle = 1:3
        for iCue = 1:4
            subtightplot(3,4,iPlot)
            hold on
            switch iStyle
                case 1
                    histfit(BEHAV.Prop.(Type{iType}).(cue_type{iCue})(:,1),20)
                    title(cueT{iCue})
                    if iCue == 1
                    ylabel('Prop of first part of session')  
                    end
%                     plot([0.5 0.5],[0 25],'g')
                    xlim([0 1])
                case 2
                    histfit(BEHAV.Prop.(Type{iType}).(cue_type{iCue})(:,last(iType)),20)
                    if iCue == 1
                    ylabel('Prop of last part of session')
                    end
%                     plot([0.5 0.5],[0 25],'g')
                    xlim([0 1])
                case 3
                histfit(BEHAV.Prop.(Type{iType}).Diff.(cue_type{iCue}),20)
                if iCue == 1
                 ylabel('Difference (last - first)')
                end
                 plot([0 0],[0 25],'g')
                 xlim([-.5 .5])
            end
            iPlot = iPlot + 1;
            box off                       
            ylim([0 25])
                  set(gca,'XTick',[],'YTick',[],'FontSize',16);
        end
    end
end


%%
for iType = 1:2
figure
block_order = [1 3 2 4];
for iPlot = 1:4
    subtightplot(2,2,iPlot)
     histfit(BEHAV.Corr.Prop.(Type{iType})(:,block_order(iPlot)),20)
     box off
     ylim([0 10])
                  set(gca,'FontSize',16);
                  if iPlot == 1
                      ylabel('1st half of block')
                      title('Block 1')
                      set(gca,'XTick',[]);
                  elseif iPlot == 3
                      ylabel('2nd half of block')
                      xlabel('Proportion correct')
                  elseif iPlot == 2
                      title('Block 2')
                      set(gca,'XTick',[],'YTick',[]);
                  elseif iPlot == 4
                      set(gca,'YTick',[]);
                        xlabel('Proportion correct')
                  end
     xlim([.5 1])
     
end
end

%%
sesh_num = 1:length(BEHAV.Corr.Prop.Half);
sesh_num = cat(2,sesh_num,sesh_num,sesh_num,sesh_num)';
num1(1:length(BEHAV.Corr.Prop.Half),1) = 1;
num2(1:length(BEHAV.Corr.Prop.Half),1) = 2;
block_num = cat(1,num1,num1,num2,num2);
half_num = cat(1,num1,num2,num1,num2);
data_cat = cat(1,BEHAV.Corr.Prop.Third(:,1),BEHAV.Corr.Prop.Third(:,2),BEHAV.Corr.Prop.Third(:,3),BEHAV.Corr.Prop.Third(:,4));

    Beh.tbl = table(sesh_num,block_num,half_num,data_cat,'VariableNames',{'SeshID','BlockID','HalfID','PropCorr'});
    Beh.tbl.BlockID = nominal(Beh.tbl.BlockID);
    Beh.tbl.HalfID = nominal(Beh.tbl.HalfID);
    Beh.lme = fitlme(Beh.tbl,'PropCorr~BlockID+HalfID+(1|SeshID)');
    Beh.lme_reduced = fitlme(Beh.tbl,'PropCorr~BlockID+(1|SeshID)');
    
    Beh.comparison = compare(Beh.lme_reduced,Beh.lme);