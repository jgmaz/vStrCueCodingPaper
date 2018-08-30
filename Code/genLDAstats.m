function Class_accuracy = genLDAstats(destination)
% function Class_accuracy = genLDAstats(destination)
%
%
% INPUTS:
%
% OUTPUTS:

cd(destination)
load(strcat(destination,'LDA_DATA_performance.mat'));
Class_accuracy.DATA = Class_accuracy_DATA;

load(strcat(destination,'LDA_SHUFF_performance.mat'));
Class_accuracy.SHUFF = Class_accuracy_SHUFF;

%% Z-scores away from shuffle (value - shuff mean / shuff std)
mdl_identifier = {'Modality','Location','Outcome'};
for iMdl = 1:length(mdl_identifier);
    Class_accuracy.Table.(mdl_identifier{iMdl}).Data = mean(Class_accuracy.DATA.(mdl_identifier{iMdl}));
    Class_accuracy.Table.(mdl_identifier{iMdl}).ShuffMEAN = mean(Class_accuracy.SHUFF.(mdl_identifier{iMdl}));
    Class_accuracy.Table.(mdl_identifier{iMdl}).ShuffSTD = std(Class_accuracy.SHUFF.(mdl_identifier{iMdl}));
    for iTime = 1:length(Class_accuracy.Table.(mdl_identifier{iMdl}).Data)
        Class_accuracy.Table.(mdl_identifier{iMdl}).Zscore(iTime) = (Class_accuracy.Table.(mdl_identifier{iMdl}).Data(iTime) - Class_accuracy.Table.(mdl_identifier{iMdl}).ShuffMEAN(iTime)) / Class_accuracy.Table.(mdl_identifier{iMdl}).ShuffSTD(iTime);
        if Class_accuracy.Table.(mdl_identifier{iMdl}).Zscore(iTime) > 1.96
            Class_accuracy.Table.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 1;
        else
            Class_accuracy.Table.(mdl_identifier{iMdl}).Zscore_recode(iTime) = 0;
        end
    end
end

%%
for iMdl = 1:length(mdl_identifier);
    Class_accuracy.DATA.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy.DATA.(mdl_identifier{iMdl}));
    Class_accuracy.SHUFF.MEAN.(mdl_identifier{iMdl}) = mean(Class_accuracy.SHUFF.(mdl_identifier{iMdl}));
    Class_accuracy.SHUFF.ZSTD.(mdl_identifier{iMdl}) = std(Class_accuracy.SHUFF.(mdl_identifier{iMdl})) *1.96;
end

save(cat(2,destination,'LDA_plotting_stats.mat'),'Class_accuracy')

end