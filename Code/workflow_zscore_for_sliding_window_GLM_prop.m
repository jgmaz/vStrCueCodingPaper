%% Z-scores away from shuffle (value - shuff mean / shuff std)

Epochs = {'cueon' 'NP' 'outcome'};
for iEpoch = 1:length(Epochs);
    if iEpoch == 1
        Predictors = {'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous'};
    else
        Predictors = {'Modality' 'Location' 'Outcome'};
    end
    
for iPred = 1:length(Predictors);
    Table.(Epochs{iEpoch}).(Predictors{iPred}).Data = GLM_window.(Epochs{iEpoch}).prop(iPred+1,:);
    Table.(Epochs{iEpoch}).(Predictors{iPred}).ShuffMEAN = mean(GLM_window_SHUFF.(Epochs{iEpoch}).prop.ALL.(Predictors{iPred}));
    Table.(Epochs{iEpoch}).(Predictors{iPred}).ShuffSTD = std(GLM_window_SHUFF.(Epochs{iEpoch}).prop.ALL.(Predictors{iPred}));
    for iTime = 1:length(Table.(Epochs{iEpoch}).(Predictors{iPred}).Data)
         Table.(Epochs{iEpoch}).(Predictors{iPred}).Zscore(iTime) = (Table.(Epochs{iEpoch}).(Predictors{iPred}).Data(iTime) - Table.(Epochs{iEpoch}).(Predictors{iPred}).ShuffMEAN(iTime)) / Table.(Epochs{iEpoch}).(Predictors{iPred}).ShuffSTD(iTime);
end
end
end

%%
figure
colors = {'b' 'r' 'y'};
for iEpoch = 1:length(Epochs)
    subplot(1,3,iEpoch)
    hold on
for iPred = 1:length(Predictors);
    plot(-.5:.1:.5,Table.(Epochs{iEpoch}).(Predictors{iPred}).Zscore,'color',colors{iPred})
end
title(Epochs{iEpoch})
plot(-.5:.05:.5,1.96,'.','color','black');
ylabel('Zscore')
xlabel('Time')
end