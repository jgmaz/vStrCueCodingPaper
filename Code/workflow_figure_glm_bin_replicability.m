% Different windows to run through: 500 bin & 100 step, 500 bin & 50 step,
% 250 bin & 50 step, 200 bin & 25 step, 
% NP and outcome has all but 500 bin & 50 step, 
% cue off only has doesn't have 200 bin & 25 step

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
Round = {'R1' 'R2' 'R3' 'R4' 'R5'};
for iEpoch = 1%:length(Epoch)
disp(iEpoch)
switch iEpoch
    case 1
mat_files = dir('2018-03-24-GLM_cueon*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
Analysis = {'.mat' '-round1.mat' '-round2.mat' '-round3.mat' '-round4.mat'};
    case 2
mat_files = dir('2018-03-26-GLM_NP*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
Analysis = {'-round1.mat' '-round2.mat' '-round3.mat' '-round4.mat' '-round5.mat'};
    case 3
mat_files = dir('2018-03-26-GLM_outcome*');
predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
Analysis = {'-round1.mat' '-round2.mat' '-round3.mat' '-round4.mat' '-round5.mat'};
end
% for iGLM = 1:length(file_order)
for iRound = 1%:5
    GLM_count = 1;
    disp(cat(2,'Round ',num2str(iRound)))
for iWindow = -.5:.1:.5
    disp(cat(2,'Window ',num2str(iWindow)))   
    if iEpoch == 1
   current_file = cat(2,'2018-03-24-GLM_',Epoch{iEpoch},'_',num2str(iWindow),Analysis{iRound});     
    else
    current_file = cat(2,'2018-03-26-GLM_',Epoch{iEpoch},'_',num2str(iWindow),Analysis{iRound});
       
    end
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
                GLM_window.(Epoch{iEpoch}).(Round{iRound}).window(iList,GLM_count) = summary_var.All.(predictor_list{iList});
                if iList > 1 && iList < 12
                    GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(iList-1,GLM_count) = nanmean(temp_Rsquared(:,iList-1));
                    GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(iList-1,GLM_count) = nanstd(temp_Rsquared(:,iList-1))/sqrt(numel(temp_Rsquared(:,iList-1))-sum(isnan(temp_Rsquared(:,iList-1))));
                end
            end
            GLM_count = GLM_count + 1;
        end
    end
end

GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop = GLM_window.(Epoch{iEpoch}).(Round{iRound}).window / 133;

end
end
%%
for iEpoch = 2:3%1:length(Epoch)
    for iRound = 1:2%:5
figure(iEpoch);
subplot(2,5,iRound)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(2:4,:))
hold on
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).prop.MEAN.Modality,GLM_window_SHUFF.(Epoch{iEpoch}).prop.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).prop.MEAN.Location,GLM_window_SHUFF.(Epoch{iEpoch}).prop.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).prop.MEAN.Outcome,GLM_window_SHUFF.(Epoch{iEpoch}).prop.SEM.Outcome,'--y',1);
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel(cat(2,'Time relative to ',Epoch{iEpoch}));

subplot(2,5,iRound+5)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(3,:),'-y',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.MEAN.Modality,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.MEAN.Location,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.25:.1:.75,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.MEAN.Outcome,GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.SEM.Outcome,'--y',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel(cat(2,'Time relative to ',Epoch{iEpoch}));

if iEpoch == 1
figure(4);
subplot(2,5,iRound)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(5:6,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('Behavioral measures'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time relative to cue-onset');

subplot(2,5,iRound+5)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(4,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(5,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(5,:),'-r',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset'); 

figure(5);
subplot(2,5,iRound)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(7:8,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
title('Task history'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time relative to cue-onset');

subplot(2,5,iRound+5)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(6,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(7,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(7,:),'-r',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');    
end

    end
end
%% cue on
figure;
subplot(2,3,1)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(2:4,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
legend({'identity' 'location' 'outcome'}); title('Cue features'); ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('Time relative to cue-onset');
subplot(2,3,2)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(5:6,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Behavioral measures'); xlabel('Time relative to cue-onset');
legend({'approach' 'latency'});
subplot(2,3,3)
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(7:8,:))
hold on
plot(.2,0.01:.01:.5,'.k'); plot(-.2,0.01:.01:.5,'.k');
ylim([0 .5]); title('Task history'); xlabel('Time relative to cue-onset');
legend({'trial number' 'previous trial'});
subplot(2,3,4)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(1,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(3,:),'-y',1);
ylim([1 8]); ylabel('Mean variance explained'); xlabel('Time relative to cue-onset');
subplot(2,3,5)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(4,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(4,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(5,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(5,:),'-r',1);
ylim([1 8]); xlabel('Time relative to cue-onset');
subplot(2,3,6)
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(6,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(6,:),'-b',1);
hold on
plot(.2,1:.1:8,'.k'); plot(-.2,1:.1:8,'.k');
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(7,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(7,:),'-r',1);
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
plot(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).prop(2:4,:))
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
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(1,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(1,:),'-b',1);
hold on;
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(2,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.25:.1:.75,GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.MEAN(3,:),GLM_window.(Epoch{iEpoch}).(Round{iRound}).Rsquared.SEM(3,:),'-y',1);
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