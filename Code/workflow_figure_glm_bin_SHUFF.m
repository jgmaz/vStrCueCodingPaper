num_Shuffs = 100;
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
for iEpoch = 2%:length(Epoch)
disp(iEpoch)
switch iEpoch
    case 1
mat_files = dir('2018-03-17-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};% ...
  %  'ModxLoc' 'ModxOut' 'LocxOut' 'OutxApp' 'ModxLocxOut'};
file_order = {'2018-03-17-GLM_cueon_-0.5-SHUFF.mat' '2018-03-17-GLM_cueon_-0.4-SHUFF.mat' ...
    '2018-03-17-GLM_cueon_-0.3-SHUFF.mat' '2018-03-17-GLM_cueon_-0.2-SHUFF.mat' ...
    '2018-03-17-GLM_cueon_-0.1-SHUFF.mat' '2018-03-17-GLM_cueon_0-SHUFF.mat' ...
    '2018-03-17-GLM_cueon_0.1-SHUFF.mat' '2018-03-17-GLM_cueon_0.2-SHUFF.mat' ...
    '2018-03-17-GLM_cueon_0.3-SHUFF.mat' '2018-03-17-GLM_cueon_0.4-SHUFF.mat' '2018-03-17-GLM_cueon_0.5-SHUFF.mat'};
    case 2
mat_files = dir('2018-03-17-GLM_NP*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-17-GLM_NP_-0.5-SHUFF.mat' '2018-03-17-GLM_NP_-0.4-SHUFF.mat' '2018-03-17-GLM_NP_-0.3-SHUFF.mat' ...
    '2018-03-17-GLM_NP_-0.2-SHUFF.mat' '2018-03-17-GLM_NP_-0.1-SHUFF.mat' '2018-03-17-GLM_NP_0-SHUFF.mat' '2018-03-17-GLM_NP_0.1-SHUFF.mat' ...
    '2018-03-17-GLM_NP_0.2-SHUFF.mat' '2018-03-17-GLM_NP_0.3-SHUFF.mat' '2018-03-17-GLM_NP_0.4-SHUFF.mat' '2018-03-17-GLM_NP_0.5-SHUFF.mat'};
    case 3
mat_files = dir('2018-03-17-GLM_outcome*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-17-GLM_outcome_-0.5-SHUFF.mat' '2018-03-17-GLM_outcome_-0.4-SHUFF.mat' '2018-03-17-GLM_outcome_-0.3-SHUFF.mat' ...
    '2018-03-17-GLM_outcome_-0.2-SHUFF.mat' '2018-03-17-GLM_outcome_-0.1-SHUFF.mat' '2018-03-17-GLM_outcome_0-SHUFF.mat' '2018-03-17-GLM_outcome_0.1-SHUFF.mat' ...
    '2018-03-17-GLM_outcome_0.2-SHUFF.mat' '2018-03-17-GLM_outcome_0.3-SHUFF.mat' '2018-03-17-GLM_outcome_0.4-SHUFF.mat' '2018-03-17-GLM_outcome_0.5-SHUFF.mat'};
    case 4
mat_files = dir('2018-03-17-GLM_cueoff*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-17-GLM_cueoff_-0.5-SHUFF.mat' '2018-03-17-GLM_cueoff_-0.4-SHUFF.mat' '2018-03-17-GLM_cueoff_-0.3-SHUFF.mat' ...
    '2018-03-17-GLM_cueoff_-0.2-SHUFF.mat' '2018-03-17-GLM_cueoff_-0.1-SHUFF.mat' '2018-03-17-GLM_cueoff_0-SHUFF.mat' '2018-03-17-GLM_cueoff_0.1-SHUFF.mat' ...
    '2018-03-17-GLM_cueoff_0.2-SHUFF.mat' '2018-03-17-GLM_cueoff_0.3-SHUFF.mat' '2018-03-17-GLM_cueoff_0.4-SHUFF.mat' '2018-03-17-GLM_cueoff_0.5-SHUFF.mat'};
end
for iGLM = 1:length(file_order)
    for iFind = 1:length(mat_files)
        if strcmp(mat_files(iFind).name,file_order{iGLM}) == 1
            load(mat_files(iFind).name,'GLM_matrices','summary_var');
            disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files))));
            for iShuff = 1:num_Shuffs
            temp_Rsquared = GLM_matrices.(cat(2,'shuff_',num2str(iShuff))).Rsquared.ALL;
            temp_Rsquared(temp_Rsquared == 0) = NaN;
            temp_Rsquared(temp_Rsquared > 50) = NaN;
            temp_Rsquared(temp_Rsquared < -50) = NaN;
            temp_Rsquared = abs(temp_Rsquared);
            for iList = 1:length(predictor_list)
                GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).window(iList,iGLM) = summary_var.(cat(2,'shuff_',num2str(iShuff))).All.(predictor_list{iList});
                if iList > 1 && iList < 12
                    GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.MEAN(iList-1,iGLM) = nanmean(temp_Rsquared(:,iList-1));
                    GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.SEM(iList-1,iGLM) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                end
            end
            end
        end
    end
end


end
%%
for iEpoch = 2%:3
for iShuff = 1:num_Shuffs
    GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).prop = GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).window / 133;
    for iList = 2:length(predictor_list)
GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})(iShuff,:) = GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).prop(iList,:);
GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList})(iShuff,:) = GLM_window.(cat(2,'shuff_',num2str(iShuff))).(Epoch{iEpoch}).Rsquared.MEAN(iList-1,:);    
    end
end
for iList = 2:length(predictor_list)
GLM_window.SHUFF.(Epoch{iEpoch}).prop.MEAN.(predictor_list{iList}) = nanmean(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList}));
GLM_window.SHUFF.(Epoch{iEpoch}).prop.SEM.(predictor_list{iList}) = nanstd(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList}))/ ...
    sqrt(numel(GLM_window.SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})(:,1))); %didn't remove nans here as i have in previous nanstd uses
GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.MEAN.(predictor_list{iList}) = nanmean(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList}));
GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.SEM.(predictor_list{iList}) = nanstd(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList}))/ ...
    sqrt(numel(GLM_window.SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList})(:,1))); %didn't remove nans here as i have in previous nanstd uses
end
end

GLM_window_SHUFF = GLM_window.SHUFF;
%% cue on w SHUFF
figure;
subplot(2,3,1)
plot(-.25:.1:.75,GLM_window.cueon.prop(2:4,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,2)
plot(-.25:.1:.75,GLM_window.cueon.prop(5:6,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Approach,GLM_window_SHUFF.cueon.prop.SEM.Approach,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Latency,GLM_window_SHUFF.cueon.prop.SEM.Latency,'--r',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:.1:.75,GLM_window.cueon.prop(7:8,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Trial,GLM_window_SHUFF.cueon.prop.SEM.Trial,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Previous,GLM_window_SHUFF.cueon.prop.SEM.Previous,'--r',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window_SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Location,GLM_window_SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window_SHUFF.cueon.Rsquared.SEM.Outcome,'--y',1);
ylim([1 8]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,5)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(4,:),GLM_window.cueon.Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(5,:),GLM_window.cueon.Rsquared.SEM(5,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Approach,GLM_window_SHUFF.cueon.Rsquared.SEM.Approach,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Latency,GLM_window_SHUFF.cueon.Rsquared.SEM.Latency,'--r',1);
ylim([1 8]); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,6)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(6,:),GLM_window.cueon.Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(7,:),GLM_window.cueon.Rsquared.SEM(7,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Trial,GLM_window_SHUFF.cueon.Rsquared.SEM.Trial,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Previous,GLM_window_SHUFF.cueon.Rsquared.SEM.Previous,'--r',1);
ylim([1 8]); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);

%% all epochs w SHUFF
figure;
subplot(2,3,1)
plot(-.25:.1:.75,GLM_window.cueon.prop(2:4,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('Cue features aligned to cue-onset'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,2)
plot(-.25:.1:.75,GLM_window.NP.prop(2:4,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.prop.MEAN.Modality,GLM_window_SHUFF.NP.prop.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.prop.MEAN.Location,GLM_window_SHUFF.NP.prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.prop.MEAN.Outcome,GLM_window_SHUFF.NP.prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
title('Cue features aligned to nosepoke'); % ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time from nosepoke (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,3)
plot(-.25:.1:.75,GLM_window.outcome.prop(2:4,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.prop.MEAN.Modality,GLM_window_SHUFF.outcome.prop.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.prop.MEAN.Location,GLM_window_SHUFF.outcome.prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.prop.MEAN.Outcome,GLM_window_SHUFF.outcome.prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
title('Cue features aligned to outcome'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time from outcome (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,4)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window_SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Location,GLM_window_SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window_SHUFF.cueon.Rsquared.SEM.Outcome,'--y',1);
ylim([1 10]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,5)
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(1,:),GLM_window.NP.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(2,:),GLM_window.NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(3,:),GLM_window.NP.Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.Rsquared.MEAN.Modality,GLM_window_SHUFF.NP.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.Rsquared.MEAN.Location,GLM_window_SHUFF.NP.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.NP.Rsquared.MEAN.Outcome,GLM_window_SHUFF.NP.Rsquared.SEM.Outcome,'--y',1);
ylim([1 10]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);
subplot(2,3,6)
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(1,:),GLM_window.outcome.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(2,:),GLM_window.outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(3,:),GLM_window.outcome.Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.Rsquared.MEAN.Modality,GLM_window_SHUFF.outcome.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.Rsquared.MEAN.Location,GLM_window_SHUFF.outcome.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.outcome.Rsquared.MEAN.Outcome,GLM_window_SHUFF.outcome.Rsquared.SEM.Outcome,'--y',1);
ylim([1 10]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)'); set(gca,'FontSize',18); xlim([-.25 .75]);

%%
for iList = 2:length(predictor_list)
GLM_window_norm.cueon.prop.MEAN.(predictor_list{iList}) = GLM_window.cueon.prop(iList,:) ./ GLM_window_SHUFF.cueon.prop.MEAN.(predictor_list{iList});
GLM_window_norm.cueon.Rsquared.MEAN.(predictor_list{iList}) = GLM_window.cueon.Rsquared.MEAN(iList,:) ./ GLM_window_SHUFF.cueon.Rsquared.MEAN.(predictor_list{iList});
GLM_window_norm.cueon.Rsquared.SEM.(predictor_list{iList}) = GLM_window.cueon.Rsquared.SEM(iList,:) ./ GLM_window_SHUFF.cueon.Rsquared.MEAN.(predictor_list{iList});
end

%% cue onset normalized to SHUFF
figure;
subplot(2,3,1)
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Modality,'-b');
hold on;
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Location,'-r');
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Outcome,'-y');
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 4]); xlabel('Time from cue-onset (s)');
subplot(2,3,2)
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Approach,'-b');
hold on;
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Latency,'-r');
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
ylim([0 4]); title('Behavioral measures'); xlabel('Time from cue-onset (s)');
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Trial,'-b');
hold on;
plot(-.25:.1:.75,GLM_window_norm.cueon.prop.MEAN.Previous,'-r');
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
ylim([0 4]); title('Task history'); xlabel('Time from cue-onset (s)');
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Modality,GLM_window_norm.cueon.Rsquared.SEM.Modality,'-b',1);
hold on
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Location,GLM_window_norm.cueon.Rsquared.SEM.Location,'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Outcome,GLM_window_norm.cueon.Rsquared.SEM.Outcome,'-y',1);
ylim([0 4]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)');
subplot(2,3,5)
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Approach,GLM_window_norm.cueon.Rsquared.SEM.Approach,'-b',1);
hold on
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Latency,GLM_window_norm.cueon.Rsquared.SEM.Latency,'-r',1);
ylim([0 4]); xlabel('Time from cue-onset (s)');
subplot(2,3,6)
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Trial,GLM_window_norm.cueon.Rsquared.SEM.Trial,'-b',1);
hold on
plot(.2,0:.1:4,'.k'); plot(-.2,0:.1:4,'.k'); plot([-.4 .8],[1 1],'-k'); xlim([-.25 .75]);
shadedErrorBar(-.25:.1:.75,GLM_window_norm.cueon.Rsquared.MEAN.Previous,GLM_window_norm.cueon.Rsquared.SEM.Previous,'-r',1);
ylim([0 4]); xlabel('Time from cue-onset (s)');

%% cue on
figure;
subplot(2,3,1)
shadedErrorBar(-.25:.1:.75,GLM_window.SHUFF.cueon.prop.MEAN.Modality,GLM_window.SHUFF.cueon.prop.SEM.Modality,'--b',1);
hold on
shadedErrorBar(-.25:.1:.75,GLM_window.SHUFF.cueon.prop.MEAN.Location,GLM_window.SHUFF.cueon.prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.SHUFF.cueon.prop.MEAN.Outcome,GLM_window.SHUFF.cueon.prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time from cue-onset (s)');
subplot(2,3,2)
plot(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.prop(5:6,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('Time from cue-onset (s)');
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.prop(7:8,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('Time from cue-onset (s)');
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(1,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(2,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(3,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(3,:),'-y',1);
ylim([1 8]); ylabel('Percent improvement to R-Squared'); xlabel('Time from cue-onset (s)');
subplot(2,3,5)
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(4,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(5,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(5,:),'-r',1);
ylim([1 8]); xlabel('Time from cue-onset (s)');
subplot(2,3,6)
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(6,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:10,'.k');
shadedErrorBar(-.25:.05:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(7,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(7,:),'-r',1);
ylim([1 8]); xlabel('Time from cue-onset (s)');

%% cue on
figure;
subplot(2,3,1)
plot(-.25:.1:.75,sum(GLM_window.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN(1:3,:)))
subplot(2,3,2)
plot(-.25:.1:.75,sum(GLM_window.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN(5:6,:)))
subplot(2,3,3)
plot(-.25:.1:.75,sum(GLM_window.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN(7:8,:)))
subplot(2,3,4)
plot(-.25:.1:.75,sum(GLM_window.(cat(2,'shuff_',num2str(iShuff))).Rsquared.MEAN))

%% NP
figure;
subplot(2,4,1)
plot(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue-onset'); ylabel('Proportion of cue-modulated units')
ylim([.1 .5]); xlabel('Time from cue-onset (s)');
subplot(2,4,2)
plot(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('NP');
ylim([.1 .5]); xlabel('Time relative to NP');
subplot(2,4,3)
plot(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Outcome');
ylim([.1 .5]); xlabel('Time relative to outcome');
subplot(2,4,4)
plot(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue-offset');
ylim([.1 .5]); xlabel('Time relative to cue-offset');
subplot(2,4,5)
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(1,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(2,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.MEAN(3,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueon.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time from cue-onset (s)');
subplot(2,4,6)
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.MEAN(1,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.MEAN(2,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.MEAN(3,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).NP.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to NP');
subplot(2,4,7)
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.MEAN(1,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.MEAN(2,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.MEAN(3,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).outcome.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to outcome');
subplot(2,4,8)
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.MEAN(1,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.MEAN(2,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.MEAN(3,:),GLM_window.(cat(2,'shuff_',num2str(iShuff))).cueoff.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to cue-offset');