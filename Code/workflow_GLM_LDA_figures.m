%% Z-scores away from shuffle (value - shuff mean / shuff std)

Epochs = {'cueon' 'NP' 'outcome'};
for iEpoch = 1:length(Epochs);
    if iEpoch == 1
        Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
        Predictors = {'Modality' 'Location' 'Outcome'};
    end
    
    for iPred = 1:length(Predictors);
        Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Data = GLM_window.(Epochs{iEpoch}).prop(iPred+1,:);
        Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).ShuffMEAN = mean(GLM_window_SHUFF.(Epochs{iEpoch}).prop.ALL.(Predictors{iPred}));
        Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).ShuffSTD = std(GLM_window_SHUFF.(Epochs{iEpoch}).prop.ALL.(Predictors{iPred}));
        for iTime = 1:length(Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Data)
            Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Zscore(iTime) = (Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Data(iTime) - Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).ShuffMEAN(iTime)) / Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).ShuffSTD(iTime);
            if Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Zscore(iTime) > 1.96
                Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Zscore_recode(iTime) = 1;
            else
                Table.GLM.(Epochs{iEpoch}).(Predictors{iPred}).Zscore_recode(iTime) = 0;
            end
        end
    end
end

%%
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'};

for iEpoch = 1:3
    if iEpoch == 1
        predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
        predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
    end
for iList = 2:length(predictor_list)
GLM_window_SHUFF.(Epoch{iEpoch}).prop.ZSTD.(predictor_list{iList}) = nanstd(GLM_window_SHUFF.(Epoch{iEpoch}).prop.ALL.(predictor_list{iList})) *1.96;
% GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.ZSTD.(predictor_list{iList}) = nanstd(GLM_window_SHUFF.(Epoch{iEpoch}).Rsquared.ALL.(predictor_list{iList})) *1.96;
end
end

%% Z-scores away from shuffle (value - shuff mean / shuff std)
mdl_identifier = {'Modality','Location','Outcome','ModxOut','ModxLoc','LocxOut','ModxLocxOut'};
for iMdl = 1:length(mdl_identifier);
    Table.LDA.(mdl_identifier{iMdl}).Data = mean(Class_accuracy.(mdl_identifier{iMdl}));
    Table.LDA.(mdl_identifier{iMdl}).ShuffMEAN = mean(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    Table.LDA.(mdl_identifier{iMdl}).ShuffSTD = std(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    for iTime = 1:length(Table.LDA.(mdl_identifier{iMdl}).Data)
         Table.LDA.(mdl_identifier{iMdl}).Zscore(iTime) = (Table.LDA.(mdl_identifier{iMdl}).Data(iTime) - Table.LDA.(mdl_identifier{iMdl}).ShuffMEAN(iTime)) / Table.LDA.(mdl_identifier{iMdl}).ShuffSTD(iTime);
            if Table.LDA.(mdl_identifier{iMdl}).Zscore(iTime) > 1.96
                Table.LDA.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 1;
            else
                Table.LDA.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 0;
            end
    end
end

%%
for iMdl = 1:length(mdl_identifier);
    Class_accuracy.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy.(mdl_identifier{iMdl}));
%     Class_accuracy.SEM.(mdl_identifier{iMdl}) = std(Class_accuracy.(mdl_identifier{iMdl}))%/sqrt(length(Class_accuracy.(mdl_identifier{iMdl})));
    Class_accuracy_SHUFF.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy_SHUFF.(mdl_identifier{iMdl}));
    Class_accuracy_SHUFF.ZSTD.(mdl_identifier{iMdl}) = std(Class_accuracy_SHUFF.(mdl_identifier{iMdl})) *1.96;
end

%% GLM cue onset
figure;
subplot(2,3,1)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
plot(-.5:.1:.5,GLM_window.cueon.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)

for iTime = 1:length(Table.GLM.cueon.Outcome.Zscore_recode)
    if Table.GLM.cueon.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.cueon.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.cueon.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
% title('Cue features aligned to cue-onset'); 
ylabel('Proportion of units')
ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

%% LDA cue onset
% figure;
% Mdls = {[1 4 5 7] [2 5 6 7] [3 4 6 7]};
Mdls = {[1 2 3] [4 5 6] [7]};
% colors = {{'b' 'r' 'y' 'g'}}; % {'r' 'b' 'y' 'g'} {'y' 'b' 'r' 'g'}};
% colors = {'b' 'r' 'g' 'c'};
colors = {[0 0 1] [1 0 0] [.3718 .7176 .3612]};
start_pt = [.08 .08 .04];

for iMdl = 1%:length(Mdls)

subplot(2,3,2)
% hold on

for iPlot = 1:length(Mdls{iMdl})

% shadedErrorBar(-.5:.1:.5,Class_accuracy_SHUFF.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy_SHUFF.ZSTD.(mdl_identifier{Mdls{iMdl}(iPlot)}),strcat('--',colors{iPlot}),1);
shadedErrorBar(-.5:.1:.5,Class_accuracy_SHUFF.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy_SHUFF.ZSTD.(mdl_identifier{Mdls{iMdl}(iPlot)}),{'--','color',colors{iPlot}},1);

hold on
% shadedErrorBar(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy.SEM.(mdl_identifier{Mdls{iMdl}(iPlot)}),colors{iPlot},1);
% plot(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),strcat('-o',colors{iPlot}),'MarkerFaceColor',colors{iPlot},'MarkerSize',4)
plot(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),'-o','color',colors{iPlot},'MarkerFaceColor',colors{iPlot},'MarkerSize',4)

for iTime = 1:length(Table.LDA.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode)
    if Table.LDA.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[start_pt(iMdl)-(iPlot*.02) start_pt(iMdl)-(iPlot*.02)], '-k', 'LineWidth',2,'color',colors{iPlot})
    end
end
end
plot(-.05,0.01:.01:1,'.k'); plot(-.45,0.01:.01:1,'.k');
% legend({'identity' 'location' 'outcome'}); 
switch iMdl
    case 1        
% title('Cue features');
ylabel('Classification rate')
    case 2
        title('Multiple cue features')
    case 3
        title('All cue features')
end

ylim([0 1]); xlabel('LDA start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);

end
box off

%% GLM other epochs
subplot(2,3,2)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Modality,GLM_window_SHUFF.NP.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Location,GLM_window_SHUFF.NP.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Outcome,GLM_window_SHUFF.NP.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);

plot(-.5:.1:.5,GLM_window.NP.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.NP.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.NP.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
for iTime = 1:length(Table.GLM.NP.Outcome.Zscore_recode)
    if Table.GLM.NP.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.NP.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.NP.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
% title('Cue features aligned to nosepoke'); % 
ylabel('Proportion of units')
ylim([0 .5]); xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

subplot(2,3,3)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Modality,GLM_window_SHUFF.outcome.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Location,GLM_window_SHUFF.outcome.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Outcome,GLM_window_SHUFF.outcome.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);

plot(-.5:.1:.5,GLM_window.outcome.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.outcome.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.outcome.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
for iTime = 1:length(Table.GLM.outcome.Outcome.Zscore_recode)
    if Table.GLM.outcome.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.outcome.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.outcome.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
% title('Cue features aligned to outcome'); %ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

%%
figure;
% subplot(2,3,1)
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.ZSTD.Modality,'--b',1);
% 
% hold on
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.ZSTD.Location,'--r',1);
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.ZSTD.Outcome,'--g',1);
% plot(-.5:.1:.5,GLM_window.cueon.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
% plot(-.5:.1:.5,GLM_window.cueon.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
% plot(-.5:.1:.5,GLM_window.cueon.prop(4,:),'-og','MarkerFaceColor','g','MarkerSize',4)
% 
% for iTime = 1:length(Table.GLM.cueon.Outcome.Zscore_recode)
%     if Table.GLM.cueon.Modality.Zscore_recode(iTime) == 1
%         plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
%     end
%     if Table.GLM.cueon.Location.Zscore_recode(iTime) == 1
%         plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
%     end
%     if Table.GLM.cueon.Outcome.Zscore_recode(iTime) == 1
%         plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','g')
%     end
% end
% plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% % legend({'identity' 'location' 'outcome'}); 
% % title('Cue features aligned to cue-onset'); 
% ylabel('Proportion of units')
% ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
% box off

subplot(2,3,2)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Modality,GLM_window_SHUFF.NP.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Location,GLM_window_SHUFF.NP.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.prop.MEAN.Outcome,GLM_window_SHUFF.NP.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);

plot(-.5:.1:.5,GLM_window.NP.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.NP.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.NP.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
for iTime = 1:length(Table.GLM.NP.Outcome.Zscore_recode)
    if Table.GLM.NP.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.NP.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.NP.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
% title('Cue features aligned to nosepoke'); % 
ylabel('Proportion of units')
ylim([0 .5]); xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

subplot(2,3,3)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Modality,GLM_window_SHUFF.outcome.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Location,GLM_window_SHUFF.outcome.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.prop.MEAN.Outcome,GLM_window_SHUFF.outcome.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);

plot(-.5:.1:.5,GLM_window.outcome.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.outcome.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.outcome.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
for iTime = 1:length(Table.GLM.outcome.Outcome.Zscore_recode)
    if Table.GLM.outcome.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.outcome.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.outcome.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'});
% title('Cue features aligned to outcome'); %ylabel('Proportion of cue-modulated units')
ylim([0 .5]); xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

% subplot(2,3,4)
% shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
% hold on
% plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
% shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
% shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),'-g',1);
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window_SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Location,GLM_window_SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
% shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window_SHUFF.cueon.Rsquared.SEM.Outcome,'--g',1);
% ylim([1 10]); ylabel('Improvement to R-Squared (%)'); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
% box off
subplot(2,3,5)
shadedErrorBar(-.5:.1:.5,GLM_window.NP.Rsquared.MEAN(1,:),GLM_window.NP.Rsquared.SEM(1,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.NP.Rsquared.MEAN(2,:),GLM_window.NP.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window.NP.Rsquared.MEAN(3,:),GLM_window.NP.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.Rsquared.MEAN.Modality,GLM_window_SHUFF.NP.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.Rsquared.MEAN.Location,GLM_window_SHUFF.NP.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.NP.Rsquared.MEAN.Outcome,GLM_window_SHUFF.NP.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
ylim([1 10]);  xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
ylabel('Improvement to R-Squared (%)');
box off
subplot(2,3,6)
shadedErrorBar(-.5:.1:.5,GLM_window.outcome.Rsquared.MEAN(1,:),GLM_window.outcome.Rsquared.SEM(1,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.outcome.Rsquared.MEAN(2,:),GLM_window.outcome.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window.outcome.Rsquared.MEAN(3,:),GLM_window.outcome.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.Rsquared.MEAN.Modality,GLM_window_SHUFF.outcome.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.Rsquared.MEAN.Location,GLM_window_SHUFF.outcome.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.outcome.Rsquared.MEAN.Outcome,GLM_window_SHUFF.outcome.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
ylim([1 10]);  xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

%% cue on w SHUFF
figure;
subplot(2,3,1)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Modality,GLM_window_SHUFF.cueon.prop.ZSTD.Modality,'--b',1);

hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Location,GLM_window_SHUFF.cueon.prop.ZSTD.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Outcome,GLM_window_SHUFF.cueon.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
plot(-.5:.1:.5,GLM_window.cueon.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)

for iTime = 1:length(Table.GLM.cueon.Outcome.Zscore_recode)
    if Table.GLM.cueon.Modality.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.cueon.Location.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
    end
    if Table.GLM.cueon.Outcome.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
% legend({'identity' 'location' 'outcome'}); 
% title('Cue features aligned to cue-onset'); 
ylabel('Proportion of units')
ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

subplot(2,3,2)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Approach,GLM_window_SHUFF.cueon.prop.ZSTD.Approach,'--b',1);
hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Latency,GLM_window_SHUFF.cueon.prop.ZSTD.Latency,'--r',1);
plot(-.5:.1:.5,GLM_window.cueon.prop(5,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(6,:),'-or','MarkerFaceColor','r','MarkerSize',4)
for iTime = 1:length(Table.GLM.cueon.Outcome.Zscore_recode)
    if Table.GLM.cueon.Approach.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.cueon.Latency.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
ylim([0 .5]); %title('Behavioral measures'); 
xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
% legend({'approach' 'latency'});
box off

subplot(2,3,3)
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Trial,GLM_window_SHUFF.cueon.prop.ZSTD.Trial,'--b',1);
hold on
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.prop.MEAN.Previous,GLM_window_SHUFF.cueon.prop.ZSTD.Previous,'--r',1);
plot(-.5:.1:.5,GLM_window.cueon.prop(7,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
plot(-.5:.1:.5,GLM_window.cueon.prop(8,:),'-or','MarkerFaceColor','r','MarkerSize',4)
for iTime = 1:length(Table.GLM.cueon.Outcome.Zscore_recode)
    if Table.GLM.cueon.Trial.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
    end
    if Table.GLM.cueon.Previous.Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
    end
end
plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
ylim([0 .5]); %title('Task history');
xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
% legend({'trial number' 'previous trial'}); 
box off

subplot(2,3,4)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(1,:),GLM_window.cueon.Rsquared.SEM(1,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(2,:),GLM_window.cueon.Rsquared.SEM(2,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(3,:),GLM_window.cueon.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window_SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Location,GLM_window_SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window_SHUFF.cueon.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
ylim([1 8]); ylabel('Improvement to R-Squared (%)'); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

subplot(2,3,5)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(4,:),GLM_window.cueon.Rsquared.SEM(4,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(5,:),GLM_window.cueon.Rsquared.SEM(5,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Approach,GLM_window_SHUFF.cueon.Rsquared.SEM.Approach,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Latency,GLM_window_SHUFF.cueon.Rsquared.SEM.Latency,'--r',1);
ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off
subplot(2,3,6)
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(6,:),GLM_window.cueon.Rsquared.SEM(6,:),'-b',1);
hold on
plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
shadedErrorBar(-.5:.1:.5,GLM_window.cueon.Rsquared.MEAN(7,:),GLM_window.cueon.Rsquared.SEM(7,:),'-r',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Trial,GLM_window_SHUFF.cueon.Rsquared.SEM.Trial,'--b',1);
shadedErrorBar(-.5:.1:.5,GLM_window_SHUFF.cueon.Rsquared.MEAN.Previous,GLM_window_SHUFF.cueon.Rsquared.SEM.Previous,'--r',1);
ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%% cue on w SHUFF
figure;
% Mdls = {[1 4 5 7] [2 5 6 7] [3 4 6 7]};
Mdls = {[1 2 3] [4 5 6] [7]};
% colors = {{'b' 'r' 'y' 'g'}}; % {'r' 'b' 'y' 'g'} {'y' 'b' 'r' 'g'}};
colors = {'b' 'r' 'g' 'c'};
start_pt = [.08 .08 .04];

for iMdl = 2:length(Mdls)

subplot(2,3,iMdl)
% hold on

for iPlot = 1:length(Mdls{iMdl})

shadedErrorBar(-.5:.1:.5,Class_accuracy_SHUFF.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy_SHUFF.ZSTD.(mdl_identifier{Mdls{iMdl}(iPlot)}),strcat('--',colors{iPlot}),1);
hold on
% shadedErrorBar(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),Class_accuracy.SEM.(mdl_identifier{Mdls{iMdl}(iPlot)}),colors{iPlot},1);
plot(-.5:.1:.5,Class_accuracy.MEAN.(mdl_identifier{Mdls{iMdl}(iPlot)}),strcat('-o',colors{iPlot}),'MarkerFaceColor',colors{iPlot},'MarkerSize',4)

for iTime = 1:length(Table.LDA.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode)
    if Table.LDA.(mdl_identifier{Mdls{iMdl}(iPlot)}).Zscore_recode(iTime) == 1
        plot([-.65+(iTime*.1) -.55+(iTime*.1)],[start_pt(iMdl)-(iPlot*.02) start_pt(iMdl)-(iPlot*.02)], '-k', 'LineWidth',2,'color',colors{iPlot})
    end
end
end
plot(-.05,0.01:.01:1,'.k'); plot(-.45,0.01:.01:1,'.k');
% legend({'identity' 'location' 'outcome'}); 
switch iMdl
    case 1        
title('Cue features');
ylabel('Classification rate')
    case 2
%         title('Multiple cue features')
ylabel('Classification rate')
    case 3
%         title('All cue features')
end

ylim([0 1]); xlabel('LDA start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off
end