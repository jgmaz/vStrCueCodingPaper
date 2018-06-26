Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
for iEpoch = 1%:length(Epoch)
disp(iEpoch)
switch iEpoch
    case 1
mat_files = dir('2018-05-21-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-05-21-GLM_cueon_APP_only_-0.5-round1.mat' '2018-05-21-GLM_cueon_APP_only_-0.4-round1.mat' '2018-05-21-GLM_cueon_APP_only_-0.3-round1.mat' ...
    '2018-05-21-GLM_cueon_APP_only_-0.2-round1.mat' '2018-05-21-GLM_cueon_APP_only_-0.1-round1.mat' '2018-05-21-GLM_cueon_APP_only_0-round1.mat' ...
    '2018-05-21-GLM_cueon_APP_only_0.1-round1.mat' '2018-05-21-GLM_cueon_APP_only_0.2-round1.mat' '2018-05-21-GLM_cueon_APP_only_0.3-round1.mat' ...
    '2018-05-21-GLM_cueon_APP_only_0.4-round1.mat' '2018-05-21-GLM_cueon_APP_only_0.5-round1.mat'};
    case 2
mat_files = dir('2018-03-*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-12-GLM_NP_-0.5-500bin.mat' '2018-03-12-GLM_NP_-0.4-500bin.mat' '2018-03-12-GLM_NP_-0.3-500bin.mat' ...
    '2018-03-12-GLM_NP_-0.2-500bin.mat' '2018-03-12-GLM_NP_-0.1-500bin.mat' '2018-03-12-GLM_NP_0-500bin.mat' '2018-03-12-GLM_NP_0.1-500bin.mat' ...
    '2018-03-12-GLM_NP_0.2-500bin.mat' '2018-03-12-GLM_NP_0.3-500bin.mat' '2018-03-12-GLM_NP_0.4-500bin.mat' '2018-03-12-GLM_NP_0.5-500bin.mat'};
    case 3
mat_files = dir('2018-03-*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-23-GLM_outcome_-0.5.mat' '2018-03-23-GLM_outcome_-0.4.mat' '2018-03-23-GLM_outcome_-0.3.mat' ...
    '2018-03-23-GLM_outcome_-0.2.mat' '2018-03-23-GLM_outcome_-0.1.mat' '2018-03-23-GLM_outcome_0.mat' '2018-03-23-GLM_outcome_0.1.mat' ...
    '2018-03-23-GLM_outcome_0.2.mat' '2018-03-23-GLM_outcome_0.3.mat' '2018-03-23-GLM_outcome_0.4.mat' '2018-03-23-GLM_outcome_0.5.mat'};
    case 4
mat_files = dir('2018-03-12-GLM_cueoff*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
file_order = {'2018-03-12-GLM_cueoff_-0.5.mat' '2018-03-12-GLM_cueoff_-0.4.mat' '2018-03-12-GLM_cueoff_-0.3.mat' ...
    '2018-03-12-GLM_cueoff_-0.2.mat' '2018-03-12-GLM_cueoff_-0.1.mat' '2018-03-12-GLM_cueoff_0.mat' '2018-03-12-GLM_cueoff_0.1.mat' ...
    '2018-03-12-GLM_cueoff_0.2.mat' '2018-03-12-GLM_cueoff_0.3.mat' '2018-03-12-GLM_cueoff_0.4.mat' '2018-03-12-GLM_cueoff_0.5.mat'};
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
ylim([0 .5]); xlabel('Time relative to cue-onset');
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
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');
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
legend({'identity' 'location' 'outcome'}); title('Cue-onset'); ylabel('Proportion of cue-modulated units')
ylim([.1 .5]); xlabel('Time relative to cue-onset');
subplot(2,4,2)
plot(-.25:.1:.75,GLM_window.NP.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('NP');
ylim([.1 .5]); xlabel('Time relative to NP');
% xlim([-.3 .3]);
subplot(2,4,3)
plot(-.25:.1:.75,GLM_window.outcome.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Outcome');
ylim([.1 .5]); xlabel('Time relative to outcome');
% xlim([-.3 .3]);
% subplot(2,4,4)
% plot(-.25:.1:.75,GLM_window.cueoff.prop(2:4,:))
% hold on
% plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); title('Cue-offset');
% ylim([.1 .5]); xlabel('Time relative to cue-offset');
% xlim([-.3 .3]);
subplot(2,4,5)
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to cue-onset');
subplot(2,4,6)
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(1,:),GLM_window.NP.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(2,:),GLM_window.NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(3,:),GLM_window.NP.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to NP');
% xlim([-.3 .3]);
subplot(2,4,7)
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(1,:),GLM_window.outcome.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(2,:),GLM_window.outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(3,:),GLM_window.outcome.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to outcome');
% xlim([-.3 .3]);
% subplot(2,4,8)
% shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(1,:),GLM_window.cueoff.Rsquared.SEM(1,:),'-b',1);
% hold on;
% shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(2,:),GLM_window.cueoff.Rsquared.SEM(2,:),'-r',1);
% shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(3,:),GLM_window.cueoff.Rsquared.SEM(3,:),'-y',1);
% plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
% ylim([2 10]); xlabel('Time relative to cue-offset');
% xlim([-.3 .3]);