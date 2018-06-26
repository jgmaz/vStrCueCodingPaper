Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
for iEpoch = 1:length(Epoch)
disp(iEpoch)
switch iEpoch
    case 1
%         mat_files = dir('2018-05-26-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-05-26-GLM_cueon_-0.5-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_-0.4-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_-0.3-SYNTH_no_prev_trial.mat' ...
%     '2018-05-26-GLM_cueon_-0.2-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_-0.1-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_0-SYNTH_no_prev_trial.mat' ...
%     '2018-05-26-GLM_cueon_0.1-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_0.2-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_0.3-SYNTH_no_prev_trial.mat' ...
%     '2018-05-26-GLM_cueon_0.4-SYNTH_no_prev_trial.mat' '2018-05-26-GLM_cueon_0.5-SYNTH_no_prev_trial.mat'};
%   
    mat_files = dir('2018-05-27-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-27-GLM_cueon_-0.5-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.4-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.3-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_-0.2-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.1-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_0.1-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.2-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.3-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_0.4-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.5-SYNTH_full_trial.mat'};
 

  case 2
% mat_files = dir('2018-05-25-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-05-25-GLM_cueon_-0.5-SYNTH.mat' '2018-05-25-GLM_cueon_-0.4-SYNTH.mat' '2018-05-25-GLM_cueon_-0.3-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_-0.2-SYNTH.mat' '2018-05-25-GLM_cueon_-0.1-SYNTH.mat' '2018-05-25-GLM_cueon_0-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_0.1-SYNTH.mat' '2018-05-25-GLM_cueon_0.2-SYNTH.mat' '2018-05-25-GLM_cueon_0.3-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_0.4-SYNTH.mat' '2018-05-25-GLM_cueon_0.5-SYNTH.mat'};

    mat_files = dir('2018-05-31-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-31-GLM_cueon_-0.5-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.4-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.3-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_-0.2-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.1-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_0.1-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.2-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.3-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_0.4-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.5-SYNTH_full_np.mat'};
 

        case 3
    mat_files = dir('2018-06-01-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-06-01-GLM_cueon_-0.5-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.4-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.3-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_-0.2-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.1-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_0.1-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.2-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.3-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_0.4-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.5-SYNTH_full_np_recode.mat'};
 
    case 4
%     mat_files = dir('2018-05-30-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-05-30-GLM_cueon_-0.5-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.4-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.3-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_-0.2-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.1-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_0.1-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.2-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.3-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_0.4-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.5-SYNTH_full_trial_ALL_CELLS.mat'};

    mat_files = dir('2018-06-02-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-06-02-GLM_cueon_-0.5-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.4-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.3-SYNTH_full_np_no_prev_trial.mat' ...
    '2018-06-02-GLM_cueon_-0.2-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.1-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0-SYNTH_full_np_no_prev_trial.mat' ...
    '2018-06-02-GLM_cueon_0.1-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.2-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.3-SYNTH_full_np_no_prev_trial.mat' ...
    '2018-06-02-GLM_cueon_0.4-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.5-SYNTH_full_np_no_prev_trial.mat'};
 

end
for iGLM = 1:length(file_order)
    for iFind = 1:length(mat_files)
        if strcmp(mat_files(iFind).name,file_order{iGLM}) == 1
            load(mat_files(iFind).name);
            disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files))));
            temp_Rsquared = GLM_matrices.Rsquared.ALL;
            temp_Rsquared(temp_Rsquared == 0) = NaN;
            temp_Rsquared(temp_Rsquared > 50) = NaN;
            temp_Rsquared(temp_Rsquared < -50) = NaN;
            temp_Rsquared = abs(temp_Rsquared);
            for iList = 1:length(predictor_list)
                GLM_window.(Epoch{iEpoch}).window(iList,iGLM) = summary_var.All.(predictor_list{iList});
                if iList > 1 && iList < 12
                    GLM_window.(Epoch{iEpoch}).Rsquared.MEAN(iList-1,iGLM) = nanmean(temp_Rsquared(:,iList-1));
                    GLM_window.(Epoch{iEpoch}).Rsquared.SEM(iList-1,iGLM) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                end
            end
            
        end
    end
end

GLM_window.(Epoch{iEpoch}).prop = GLM_window.(Epoch{iEpoch}).window / 133;
end
%% cue on
figure;
subplot(2,3,1)
plot(-.25:.1:.75,GLM_window.cueon.prop(2:4,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 1]); xlabel('Time relative to cue-onset');
subplot(2,3,2)
plot(-.25:.1:.75,GLM_window.cueon.prop(5:6,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('Time relative to cue-onset');
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:.1:.75,GLM_window.cueon.prop(7:8,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('Time relative to cue-onset');
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
 ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');
subplot(2,3,5)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(4,:),GLM_window.cueon.Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(5,:),GLM_window.cueon.Rsquared.SEM(5,:),'-r',1);
ylim([1 8]); xlabel('Time relative to cue-onset');
subplot(2,3,6)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(6,:),GLM_window.cueon.Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(7,:),GLM_window.cueon.Rsquared.SEM(7,:),'-r',1);
ylim([1 8]); xlabel('Time relative to cue-onset');

%% cue on
figure;
subplot(2,3,1)
plot(-.25:.1:.75,sum(GLM_window.Rsquared.MEAN(1:3,:)))
subplot(2,3,2)
plot(-.25:.1:.75,sum(GLM_window.Rsquared.MEAN(5:6,:)))
subplot(2,3,3)
plot(-.25:.1:.75,sum(GLM_window.Rsquared.MEAN(7:8,:)))
subplot(2,3,4)
plot(-.25:.1:.75,sum(GLM_window.Rsquared.MEAN))

%% NP
figure;
subplot(2,4,1)
plot(-.25:.1:.75,GLM_window.cueon.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full trial'); ylabel('Proportion of cue-modulated units')
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
subplot(2,4,2)
plot(-.25:.1:.75,GLM_window.NP.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,3)
plot(-.25:.1:.75,GLM_window.outcome.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np recode prev trial');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,4)
plot(-.25:.1:.75,GLM_window.cueoff.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np no prev trial');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,5)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 50]); xlabel('Time centered relative to cue-onset');
subplot(2,4,6)
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(1,:),GLM_window.NP.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(2,:),GLM_window.NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(3,:),GLM_window.NP.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 50]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,7)
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(1,:),GLM_window.outcome.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(2,:),GLM_window.outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(3,:),GLM_window.outcome.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 50]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,8)
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(1,:),GLM_window.cueoff.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(2,:),GLM_window.cueoff.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(3,:),GLM_window.cueoff.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 50]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);

%%
Epoch = {'d_org' 'd_recode' 'd_app' 's_best' 's_ft' 's_ft_r' 's_np' 's_np_r'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
for iEpoch = 1:length(Epoch)
disp(iEpoch)
switch iEpoch
%     case 1
%       mat_files = dir('2018-03-24-GLM_cueon*');
% file_order = {'2018-03-24-GLM_cueon_-0.5.mat' '2018-03-24-GLM_cueon_-0.4.mat' '2018-03-24-GLM_cueon_-0.3.mat' ...
%     '2018-03-24-GLM_cueon_-0.2.mat' '2018-03-24-GLM_cueon_-0.1.mat' '2018-03-24-GLM_cueon_0.mat' ...
%     '2018-03-24-GLM_cueon_0.1.mat' '2018-03-24-GLM_cueon_0.2.mat' '2018-03-24-GLM_cueon_0.3.mat' ...
%     '2018-03-24-GLM_cueon_0.4.mat' '2018-03-24-GLM_cueon_0.5.mat'};

    case 1
      mat_files = dir('2018-06-06-GLM_cueon*');
file_order = {'2018-06-06-GLM_cueon_-0.5-replication.mat' '2018-06-06-GLM_cueon_-0.4-replication.mat' '2018-06-06-GLM_cueon_-0.3-replication.mat' ...
    '2018-06-06-GLM_cueon_-0.2-replication.mat' '2018-06-06-GLM_cueon_-0.1-replication.mat' '2018-06-06-GLM_cueon_0-replication.mat' ...
    '2018-06-06-GLM_cueon_0.1-replication.mat' '2018-06-06-GLM_cueon_0.2-replication.mat' '2018-06-06-GLM_cueon_0.3-replication.mat' ...
    '2018-06-06-GLM_cueon_0.4-replication.mat' '2018-06-06-GLM_cueon_0.5-replication.mat'};
  
    case 2
          mat_files = dir('2018-06-07-GLM_cueon*');
        file_order = {'2018-06-07-GLM_cueon_-0.5-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_-0.4-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_-0.3-recodePrevTrial_v2.mat' ...
    '2018-06-07-GLM_cueon_-0.2-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_-0.1-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_0-recodePrevTrial_v2.mat' ...
    '2018-06-07-GLM_cueon_0.1-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_0.2-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_0.3-recodePrevTrial_v2.mat' ...
    '2018-06-07-GLM_cueon_0.4-recodePrevTrial_v2.mat' '2018-06-07-GLM_cueon_0.5-recodePrevTrial_v2.mat'};
  
    case 3
         mat_files = dir('2018-06-11-GLM_cueon*');
        file_order = {'2018-06-11-GLM_cueon_-0.5-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_-0.4-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_-0.3-recodePrevTrial_v2_prev_App_only.mat' ...
    '2018-06-11-GLM_cueon_-0.2-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_-0.1-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_0-recodePrevTrial_v2_prev_App_only.mat' ...
    '2018-06-11-GLM_cueon_0.1-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_0.2-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_0.3-recodePrevTrial_v2_prev_App_only.mat' ...
    '2018-06-11-GLM_cueon_0.4-recodePrevTrial_v2_prev_App_only.mat' '2018-06-11-GLM_cueon_0.5-recodePrevTrial_v2_prev_App_only.mat'};

%      mat_files = dir('2018-06-09-GLM_cueon*');
%         file_order = {'2018-06-09-GLM_cueon_-0.5-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_-0.4-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_-0.3-recodePrevBeh_Cue.mat' ...
%     '2018-06-09-GLM_cueon_-0.2-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_-0.1-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_0-recodePrevBeh_Cue.mat' ...
%     '2018-06-09-GLM_cueon_0.1-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_0.2-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_0.3-recodePrevBeh_Cue.mat' ...
%     '2018-06-09-GLM_cueon_0.4-recodePrevBeh_Cue.mat' '2018-06-09-GLM_cueon_0.5-recodePrevBeh_Cue.mat'};


%           mat_files = dir('2018-06-05-GLM_cueon*');
%         file_order = {'2018-06-05-GLM_cueon_-0.5-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_-0.4-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_-0.3-recodePrevTrial_v3.mat' ...
%     '2018-06-05-GLM_cueon_-0.2-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_-0.1-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_0-recodePrevTrial_v3.mat' ...
%     '2018-06-05-GLM_cueon_0.1-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_0.2-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_0.3-recodePrevTrial_v3.mat' ...
%     '2018-06-05-GLM_cueon_0.4-recodePrevTrial_v3.mat' '2018-06-05-GLM_cueon_0.5-recodePrevTrial_v3.mat'};

    case 4
        mat_files = dir('2018-06-04-GLM_cueon*');
        file_order = {'2018-06-04-GLM_cueon_-0.5-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_-0.4-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_-0.3-SYNTH_np_and_trial_recode.mat' ...
    '2018-06-04-GLM_cueon_-0.2-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_-0.1-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_0-SYNTH_np_and_trial_recode.mat' ...
    '2018-06-04-GLM_cueon_0.1-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_0.2-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_0.3-SYNTH_np_and_trial_recode.mat' ...
    '2018-06-04-GLM_cueon_0.4-SYNTH_np_and_trial_recode.mat' '2018-06-04-GLM_cueon_0.5-SYNTH_np_and_trial_recode.mat'};

%     case 4
% mat_files = dir('2018-05-25-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-05-25-GLM_cueon_-0.5-SYNTH.mat' '2018-05-25-GLM_cueon_-0.4-SYNTH.mat' '2018-05-25-GLM_cueon_-0.3-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_-0.2-SYNTH.mat' '2018-05-25-GLM_cueon_-0.1-SYNTH.mat' '2018-05-25-GLM_cueon_0-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_0.1-SYNTH.mat' '2018-05-25-GLM_cueon_0.2-SYNTH.mat' '2018-05-25-GLM_cueon_0.3-SYNTH.mat' ...
%     '2018-05-25-GLM_cueon_0.4-SYNTH.mat' '2018-05-25-GLM_cueon_0.5-SYNTH.mat'};

    case 5
    mat_files = dir('2018-05-28-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-28-GLM_cueon_-0.5-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_-0.4-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_-0.3-SYNTH_full_trial_no_prev_trial.mat' ...
    '2018-05-28-GLM_cueon_-0.2-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_-0.1-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_0-SYNTH_full_trial_no_prev_trial.mat' ...
    '2018-05-28-GLM_cueon_0.1-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_0.2-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_0.3-SYNTH_full_trial_no_prev_trial.mat' ...
    '2018-05-28-GLM_cueon_0.4-SYNTH_full_trial_no_prev_trial.mat' '2018-05-28-GLM_cueon_0.5-SYNTH_full_trial_no_prev_trial.mat'};
 
    
    case 6
    mat_files = dir('2018-05-27-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-27-GLM_cueon_-0.5-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.4-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.3-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_-0.2-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_-0.1-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_0.1-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.2-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.3-SYNTH_full_trial.mat' ...
    '2018-05-27-GLM_cueon_0.4-SYNTH_full_trial.mat' '2018-05-27-GLM_cueon_0.5-SYNTH_full_trial.mat'};
 

  case 7

    mat_files = dir('2018-05-31-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-31-GLM_cueon_-0.5-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.4-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.3-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_-0.2-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_-0.1-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_0.1-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.2-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.3-SYNTH_full_np.mat' ...
    '2018-05-31-GLM_cueon_0.4-SYNTH_full_np.mat' '2018-05-31-GLM_cueon_0.5-SYNTH_full_np.mat'};
 

        case 8
    mat_files = dir('2018-06-01-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-06-01-GLM_cueon_-0.5-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.4-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.3-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_-0.2-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_-0.1-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_0.1-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.2-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.3-SYNTH_full_np_recode.mat' ...
    '2018-06-01-GLM_cueon_0.4-SYNTH_full_np_recode.mat' '2018-06-01-GLM_cueon_0.5-SYNTH_full_np_recode.mat'};
 
%     mat_files = dir('2018-05-30-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-05-30-GLM_cueon_-0.5-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.4-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.3-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_-0.2-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_-0.1-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_0.1-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.2-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.3-SYNTH_full_trial_ALL_CELLS.mat' ...
%     '2018-05-30-GLM_cueon_0.4-SYNTH_full_trial_ALL_CELLS.mat' '2018-05-30-GLM_cueon_0.5-SYNTH_full_trial_ALL_CELLS.mat'};

%     mat_files = dir('2018-06-02-GLM_cueon*');
% predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
% file_order = {'2018-06-02-GLM_cueon_-0.5-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.4-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.3-SYNTH_full_np_no_prev_trial.mat' ...
%     '2018-06-02-GLM_cueon_-0.2-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_-0.1-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0-SYNTH_full_np_no_prev_trial.mat' ...
%     '2018-06-02-GLM_cueon_0.1-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.2-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.3-SYNTH_full_np_no_prev_trial.mat' ...
%     '2018-06-02-GLM_cueon_0.4-SYNTH_full_np_no_prev_trial.mat' '2018-06-02-GLM_cueon_0.5-SYNTH_full_np_no_prev_trial.mat'};
 

end
for iGLM = 1:length(file_order)
    for iFind = 1:length(mat_files)
        if strcmp(mat_files(iFind).name,file_order{iGLM}) == 1
            load(mat_files(iFind).name);
            disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files))));
            temp_Rsquared = GLM_matrices.Rsquared.ALL;
            temp_Rsquared(temp_Rsquared == 0) = NaN;
            temp_Rsquared(temp_Rsquared > 50) = NaN;
            temp_Rsquared(temp_Rsquared < -50) = NaN;
            temp_Rsquared = abs(temp_Rsquared);
            for iList = 1:length(predictor_list)
                GLM_window.(Epoch{iEpoch}).window(iList,iGLM) = summary_var.All.(predictor_list{iList});
                if iList > 1 && iList < 12
                    GLM_window.(Epoch{iEpoch}).Rsquared.MEAN(iList-1,iGLM) = nanmean(temp_Rsquared(:,iList-1));
                    GLM_window.(Epoch{iEpoch}).Rsquared.SEM(iList-1,iGLM) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                end
            end
            
        end
    end
end

GLM_window.(Epoch{iEpoch}).prop = GLM_window.(Epoch{iEpoch}).window / 133;
end

%% fig
figure;
for iEpoch = 1:length(Epoch)
subplot(2,4,iEpoch)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title(Epoch{iEpoch}); 
if iEpoch == 1 || iEpoch == 5
ylabel('Proportion of cue-modulated units')
end

ylim([0 .6]); xlabel('Time centered relative to cue-onset');
end

%%
subplot(2,4,2)
plot(-.25:.1:.75,GLM_window.NP.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,3)
plot(-.25:.1:.75,GLM_window.outcome.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np recode prev trial');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
subplot(2,4,4)
plot(-.25:.1:.75,GLM_window.cueoff.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Full np no prev trial');
ylim([0 .6]); xlabel('Time centered relative to cue-onset');
% xlim([-.3 .3]);
