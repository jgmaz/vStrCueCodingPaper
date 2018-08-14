% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\',sesh.session_id(1:4),'\Mat'));
cd('E:\Jimmie\Jimmie\Analysis\Mat');

count = 1;

mat_files = dir('*.mat');
for kk = 1:length(dir('*.mat'))
    disp(kk)
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    
%     block_drift.block1_length(kk) = length(FRATE.Cue.Trial_firing_rate_block1);
%     block_drift.block1_half(kk) = round(block_drift.block1_length(kk) / 2);
%     block_drift.b1_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)));
%     block_drift.b1_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
%     block_drift.MWU_b1(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
%     
%     block_drift.block2_length(kk) = length(FRATE.Cue.Trial_firing_rate_block2);
%     block_drift.block2_half(kk) = round(block_drift.block2_length(kk) / 2);
%     block_drift.b2_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)));
%     block_drift.b2_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
%     block_drift.MWU_b2(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
%     
%     switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
%                 case 0
    
    mat_files_drift{count} = mat_files(kk).name; 
        
    %%    
    rng('shuffle')
    num_trials_light = length(PETH.Trial.ALL.trials_light_PETH(1,:));
%     temp_idx = randperm(num_trials);
%     half_1_light = temp_idx(1:num_trials/2);
%     half_2_light = temp_idx(num_trials/2+1:end);
    
    num_trials_sound = length(PETH.Trial.ALL.trials_sound_PETH(1,:));
%     temp_idx = randperm(num_trials);
%     half_1_sound = temp_idx(1:num_trials/2);
%     half_2_sound = temp_idx(num_trials/2+1:end);
    
    for iBin = 1:15001 %generate averaged responses
%         if iBin <= length(PETH.Trial.ALL.trials_rew_PETH)
    PETHS.mod.trials_light_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_light_PETH(iBin,1:round(num_trials_light/2)));
    PETHS.mod.trials_light_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_light_PETH(iBin,round(num_trials_light/2)+1:end));
    
    PETHS.mod.trials_sound_1st_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_sound_PETH(iBin,1:round(num_trials_sound/2)));
    PETHS.mod.trials_sound_2nd_half{count}(iBin,:) = nanmean(PETH.Trial.ALL.trials_sound_PETH(iBin,round(num_trials_sound/2)+1:end));
%         end
    end
    
    switch sesh.block_order
        case 1
            PETHS.block.block1_1st_half{count} = PETHS.mod.trials_light_1st_half{count};
            PETHS.block.block1_2nd_half{count} = PETHS.mod.trials_light_2nd_half{count};
            PETHS.block.block2_1st_half{count} = PETHS.mod.trials_sound_1st_half{count};
           PETHS.block.block2_2nd_half{count} =  PETHS.mod.trials_sound_2nd_half{count};
        case 2
            PETHS.block.block1_1st_half{count} = PETHS.mod.trials_sound_1st_half{count};
            PETHS.block.block1_2nd_half{count} = PETHS.mod.trials_sound_2nd_half{count};
            PETHS.block.block2_1st_half{count} = PETHS.mod.trials_light_1st_half{count};
           PETHS.block.block2_2nd_half{count} =  PETHS.mod.trials_light_2nd_half{count};
    end
    
    temp = corrcoef(PETHS.block.block1_1st_half{count},PETHS.block.block1_2nd_half{count});
    Corr.block.block1(count) = temp(2);
    
    temp = corrcoef(PETHS.block.block2_1st_half{count},PETHS.block.block2_2nd_half{count});
    Corr.block.block2(count) = temp(2);
    
    temp = corrcoef(PETHS.block.block1_2nd_half{count},PETHS.block.block2_1st_half{count});
    Corr.block.block1v2(count) = temp(2);
            
    count = count + 1;
    
%         case 1
%     end
end

%%
[~,Corr.block.ttest.within] = ttest(Corr.block.block1,Corr.block.block2);
[~,Corr.block.ttest.across1] = ttest(Corr.block.block1,Corr.block.block1v2);
[~,Corr.block.ttest.across2] = ttest(Corr.block.block2,Corr.block.block1v2);

%%
anova_input = cat(2,Corr.block.block1',Corr.block.block2',Corr.block.block1v2');
[p_anova tbl_anova stats_anova] = anova1(anova_input);
stats_multicompare = multcompare(stats_anova);

%%
rm_table = table(Corr.block.block1',Corr.block.block1v2',Corr.block.block2',...
    'VariableNames',{'b1','b1v2','b2'});
rm_model = fitrm(rm_table,'b1-2~1');
rm_output = ranova(rm_model);

%%
figure
subplot(3,4,1)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.ctrl.sound))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v ctrl sound; p = ',num2str(Corr.mod.ttest.ctrl.light_sound)))
xlabel('ctrl light'); ylabel('ctrl sound');
subplot(3,4,2)
scatter(abs(Corr.mod.cond.l1_v_s1),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp2; p = ',num2str(Corr.mod.ttest.ctrl.cond1_2)))
xlabel('comparison 1'); ylabel('comparison 2');
subplot(3,4,3)
scatter(abs(Corr.mod.cond.l1_v_s2),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'comp2 v comp3; p = ',num2str(Corr.mod.ttest.ctrl.cond2_3)))
xlabel('comparison 2'); ylabel('comparison 3');
subplot(3,4,4)
scatter(abs(Corr.mod.cond.l1_v_s1),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'comp1 v comp4; p = ',num2str(Corr.mod.ttest.ctrl.cond1_4)))
xlabel('comparison 1'); ylabel('comparison 4');

subplot(3,4,5)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l1_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp1; p = ',num2str(Corr.mod.ttest.cond.light_cond1)))
xlabel('ctrl light'); ylabel('comparison 1');
subplot(3,4,6)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp2; p = ',num2str(Corr.mod.ttest.cond.light_cond2)))
xlabel('ctrl light'); ylabel('comparison 2');
subplot(3,4,7)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp3; p = ',num2str(Corr.mod.ttest.cond.light_cond3)))
xlabel('ctrl light'); ylabel('comparison 3');
subplot(3,4,8)
scatter(abs(Corr.mod.ctrl.light),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl light v comp4; p = ',num2str(Corr.mod.ttest.cond.light_cond4)))
xlabel('ctrl light'); ylabel('comparison 4');

subplot(3,4,9)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l1_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp1; p = ',num2str(Corr.mod.ttest.cond.sound_cond1)))
xlabel('ctrl sound'); ylabel('comparison 1');
subplot(3,4,10)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l1_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp2; p = ',num2str(Corr.mod.ttest.cond.sound_cond2)))
xlabel('ctrl sound'); ylabel('comparison 2');
subplot(3,4,11)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l2_v_s1))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp3; p = ',num2str(Corr.mod.ttest.cond.sound_cond3)))
xlabel('ctrl sound'); ylabel('comparison 3');
subplot(3,4,12)
scatter(abs(Corr.mod.ctrl.sound),abs(Corr.mod.cond.l2_v_s2))
hold on; line([0 1],[0 1]); title(cat(2,'ctrl sound v comp4; p = ',num2str(Corr.mod.ttest.cond.sound_cond4)))
xlabel('ctrl sound'); ylabel('comparison 4');

%% set up variables for LME (mod)
temp_Corr = 1:length(Corr.mod.ctrl.light);
Corr_mod_summary(:,1) = cat(2,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr,temp_Corr)';

temp_length = length(Corr.mod.ctrl.light);
for iCell = 1:temp_length
temp_Corr_cells{iCell,1} = 'A';
end
for iCell = temp_length+1:temp_length*2
temp_Corr_cells{iCell,1} = 'B';
end
for iCell = temp_length*2+1:temp_length*3
temp_Corr_cells{iCell,1} = 'C';
end
for iCell = temp_length*3+1:temp_length*4
temp_Corr_cells{iCell,1} = 'D';
end
for iCell = temp_length*4+1:temp_length*5
temp_Corr_cells{iCell,1} = 'E';
end
for iCell = temp_length*5+1:temp_length*6
temp_Corr_cells{iCell,1} = 'F';
end

Corr_mod_summary(:,2) = cat(2,Corr.mod.ctrl.light,Corr.mod.ctrl.sound,Corr.mod.cond.l1_v_s1,Corr.mod.cond.l1_v_s2,Corr.mod.cond.l2_v_s1,Corr.mod.cond.l2_v_s2)';

%% linear mixed effects model (mod)
Corr.mod.LME.tbl = table(Corr_mod_summary(:,1),temp_Corr_cells,Corr_mod_summary(:,2),'VariableNames',{'CellID','compType','corrCoeff'});
Corr.mod.LME.tbl.compType = nominal(Corr.mod.LME.tbl.compType);
Corr.mod.LME.lme = fitlme(Corr.mod.LME.tbl,'corrCoeff~compType+(1|CellID)');%,'DummyVarCoding','effects');% +(CueType-1|RatID)'); <- matlab includes this, not sure why %lme for Trial length with cue type as fixed effect and a random intercept for each rat.
Corr.mod.LME.lme_reduced = fitlme(Corr.mod.LME.tbl,'corrCoeff~1+(1|CellID)');

Corr.mod.LME.comparison = compare(Corr.mod.LME.lme_reduced,Corr.mod.LME.lme);

%% means for table
Corr.mod.mean.light = nanmean(abs(Corr.mod.ctrl.light));
Corr.mod.mean.sound = nanmean(abs(Corr.mod.ctrl.sound));
Corr.mod.mean.cond1 = nanmean(abs(Corr.mod.cond.l1_v_s1));
Corr.mod.mean.cond2 = nanmean(abs(Corr.mod.cond.l1_v_s2));
Corr.mod.mean.cond3 = nanmean(abs(Corr.mod.cond.l2_v_s1));
Corr.mod.mean.cond4 = nanmean(abs(Corr.mod.cond.l2_v_s2));

%% RM ANOVA
%table coeff, 1v2, lvs, subj)
Predictor = {'mod'};% 'loc' 'out'};
meas = {'A' 'B' 'D' 'E'}; % no C and F as they are 1v1 and 2v2, respectively.

for iPred = 1:3
    for iMeas = 1:length(meas)
        startT = find(Corr.(Predictor{iPred}).LME.tbl.compType == meas{iMeas},1,'first');
        endT = find(Corr.(Predictor{iPred}).LME.tbl.compType == meas{iMeas},1,'last');
        count = 1;
        for iData = startT:endT
            RMA.(Predictor{iPred}).RAW(count,iMeas) = Corr.(Predictor{iPred}).LME.tbl.corrCoeff(iData);
            count = count + 1;
        end
    end
    RMA.(Predictor{iPred}).tbl = table(RMA.(Predictor{iPred}).RAW(:,1),RMA.(Predictor{iPred}).RAW(:,2),RMA.(Predictor{iPred}).RAW(:,3),RMA.(Predictor{iPred}).RAW(:,4),'VariableNames',{'w1' 'w2' 'a1' 'a2'});
end

%%
f1 = [1 2 1 2];
f2 = [1 2 2 1];

for iPred = 1:3
    RMA.(Predictor{iPred}).RAWcat = [];
    RMA.(Predictor{iPred}).RAWf1 = [];
    RMA.(Predictor{iPred}).RAWf2 = [];
    RMA.(Predictor{iPred}).RAWid = [];
    for iMeas = 1:length(meas)
        temp1 = [];
        temp2 = [];
        temp3 = [];
        RMA.(Predictor{iPred}).RAWcat = cat(1,RMA.(Predictor{iPred}).RAWcat,RMA.(Predictor{iPred}).RAW(:,iMeas));
        temp1(1:length(RMA.(Predictor{iPred}).RAW(:,iMeas)),1) = f1(iMeas);
        temp2(1:length(RMA.(Predictor{iPred}).RAW(:,iMeas)),1) = f2(iMeas);
        temp3(:,1) = 1:length(RMA.(Predictor{iPred}).RAW(:,iMeas));
        RMA.(Predictor{iPred}).RAWf1 = cat(1,RMA.(Predictor{iPred}).RAWf1,temp1);
        RMA.(Predictor{iPred}).RAWf2 = cat(1,RMA.(Predictor{iPred}).RAWf2,temp2);
        RMA.(Predictor{iPred}).RAWid = cat(1,RMA.(Predictor{iPred}).RAWid,temp3);       
    end
    RMA.(Predictor{iPred}).RMtbl = table(RMA.(Predictor{iPred}).RAWcat,categorical(RMA.(Predictor{iPred}).RAWf1),categorical(RMA.(Predictor{iPred}).RAWf2),RMA.(Predictor{iPred}).RAWid,...
        'VariableNames',{'Corr' 'factor1' 'factor2' 'cellid'});
%     RMA.(Predictor{iPred}).RM = fitrm(RMA.(Predictor{iPred}).RMtbl,'Corr ~ factor1*factor2','WithinDesign',within)
end

%%
for iPred = 1:3
    [RMA.(Predictor{iPred}).RM.p RMA.(Predictor{iPred}).RM.tbl RMA.(Predictor{iPred}).RM.stats] = anovan(RMA.(Predictor{iPred}).RAWcat,{categorical(RMA.(Predictor{iPred}).RAWf1),categorical(RMA.(Predictor{iPred}).RAWf2),categorical(RMA.(Predictor{iPred}).RAWid)},...
    'model','interaction','random',[3]);
end