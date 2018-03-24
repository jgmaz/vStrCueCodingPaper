% Different windows to run through: 500 bin & 100 step, 500 bin & 50 step,
% 250 bin & 50 step, 200 bin & 25 step, 
% NP and outcome has all but 500 bin & 50 step, 
% cue off only has doesn't have 200 bin & 25 step

%%not finalized%%

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
for iEpoch = 1:length(Epoch)
disp(iEpoch)
switch iEpoch
    case 1
mat_files = dir('2018-03-12-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
Step = .100;
Analysis = {'.mat'};
    case 2
mat_files = dir('2018-03-12-GLM_NP*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
Step = .100;
Analysis = {'-500bin.mat'};
    case 3
mat_files = dir('2018-03-12-GLM_outcome*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
Step = .100;
Analysis = {'-500bin.mat'};
end
% for iGLM = 1:length(file_order)
for iStep = 1%:length(Step)
    GLM_count = 1;
    disp(cat(2,'Step ',num2str(iStep)))
for iWindow = -.5:Step(iStep):.35
    disp(cat(2,'Window ',num2str(iWindow)))
    current_file = cat(2,'2018-03-12-GLM_',Epoch{iEpoch},'_',num2str(iWindow),Analysis{iStep});
    for iFind = 1:length(mat_files)        
        if strcmp(mat_files(iFind).name,current_file) == 1
            load(mat_files(iFind).name);
            disp('file found');
            temp_Rsquared = GLM_matrices.Rsquared.ALL;
            temp_Rsquared(temp_Rsquared == 0) = NaN;
            temp_Rsquared(temp_Rsquared > 50) = NaN;
            temp_Rsquared(temp_Rsquared < -50) = NaN;
            temp_Rsquared = abs(temp_Rsquared);
            for iList = 1:length(predictor_list)
                GLM_window.(Epoch{iEpoch}).(var_name{iStep}).window(iList,GLM_count) = summary_var.All.(predictor_list{iList});
                if iList > 1 && iList < 12
                    GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(iList-1,GLM_count) = nanmean(temp_Rsquared(:,iList-1));
                    GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(iList-1,GLM_count) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                end
            end
            GLM_count = GLM_count + 1;
        end
    end
end

GLM_window.(Epoch{iEpoch}).prop = GLM_window.(Epoch{iEpoch}).(var_name{iStep}).window / 133;

end
end

%%
figure;
subplot(2,3,1)
plot(-.25:Step(iStep):.60,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(2:4,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time relative to cue-onset');
% subplot(2,3,2)
% plot(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(5:6,:))
% hold on
% plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% ylim([0 .5]); title('Behavioral measures'); xlabel('Time relative to cue-onset');
% legend({'approach' 'latency'});
% subplot(2,3,3)
% plot(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(7:8,:))
% hold on
% plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% ylim([0 .5]); title('Task history'); xlabel('Time relative to cue-onset');
% legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:Step(iStep):.60,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:Step(iStep):.60,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:Step(iStep):.60,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(3,:),'-y',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');
% subplot(2,3,5)
% shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(4,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(4,:),'-b',1);
% hold on
% plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
% shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(5,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(5,:),'-r',1);
% ylim([1 8]); xlabel('Time relative to cue-onset');
% subplot(2,3,6)
% shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(6,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(6,:),'-b',1);
% hold on
% plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
% shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(7,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(7,:),'-r',1);
% ylim([1 8]); xlabel('Time relative to cue-onset');
end
end
%% cue on
figure;
subplot(2,3,1)
plot(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(2:4,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time relative to cue-onset');
subplot(2,3,2)
plot(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(5:6,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('Time relative to cue-onset');
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(7:8,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('Time relative to cue-onset');
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(3,:),'-y',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');
subplot(2,3,5)
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(4,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(5,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(5,:),'-r',1);
ylim([1 8]); xlabel('Time relative to cue-onset');
subplot(2,3,6)
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(6,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:Step(iStep):.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(7,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(7,:),'-r',1);
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
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).prop(2:4,:))
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
subplot(2,4,3)
plot(-.25:.1:.75,GLM_window.outcome.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Outcome');
ylim([.1 .5]); xlabel('Time relative to outcome');
subplot(2,4,4)
plot(-.25:.1:.75,GLM_window.cueoff.prop(2:4,:))
hold on
plot(.2,0.1:.01:.5,'.k'); plot(-.2,0.1:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue-offset');
ylim([.1 .5]); xlabel('Time relative to cue-offset');
subplot(2,4,5)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(var_name{iStep}).Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to cue-onset');
subplot(2,4,6)
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(1,:),GLM_window.NP.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(2,:),GLM_window.NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.NP.Rsquared.MEAN(3,:),GLM_window.NP.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to NP');
subplot(2,4,7)
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(1,:),GLM_window.outcome.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(2,:),GLM_window.outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.outcome.Rsquared.MEAN(3,:),GLM_window.outcome.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to outcome');
subplot(2,4,8)
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(1,:),GLM_window.cueoff.Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(2,:),GLM_window.cueoff.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.cueoff.Rsquared.MEAN(3,:),GLM_window.cueoff.Rsquared.SEM(3,:),'-y',1);
plot(.2,2:.1:10,'.k'); plot(-.2,2:.1:10,'.k');
ylim([2 10]); xlabel('Time relative to cue-offset');