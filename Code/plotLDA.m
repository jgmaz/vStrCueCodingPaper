function LDAplot = plotLDA(Class_accuracy)
% function LDAplot = plotLDA(Class_accuracy)
%
%
% INPUTS:
%
% OUTPUTS:

%% LDA cue-onset
figure('units','normalized','outerposition',[0 0 1 1]);
mdl_identifier = {'Modality','Location','Outcome'};
colors = {[0 0 1] [1 0 0] [.3718 .7176 .3612]};
start_pt = [.08];

subplot(2,3,2)
for iPlot = 1:length(mdl_identifier)
    
    shadedErrorBar(-.5:.1:.5,Class_accuracy.SHUFF.MEAN.(mdl_identifier{iPlot}),Class_accuracy.SHUFF.ZSTD.(mdl_identifier{iPlot}),{'--','color',colors{iPlot}},1);
    
    hold on
    plot(-.5:.1:.5,Class_accuracy.DATA.MEAN.(mdl_identifier{iPlot}),'-o','color',colors{iPlot},'MarkerFaceColor',colors{iPlot},'MarkerSize',4)
    
    for iTime = 1:length(Class_accuracy.Table.(mdl_identifier{iPlot}).Zscore_recode)
        if Class_accuracy.Table.(mdl_identifier{iPlot}).Zscore_recode(iTime) == 1
            plot([-.65+(iTime*.1) -.55+(iTime*.1)],[start_pt-(iPlot*.02) start_pt-(iPlot*.02)], '-k', 'LineWidth',2,'color',colors{iPlot})
        end
    end
end
plot(-.05,0.01:.02:1,'.k'); plot(-.45,0.01:.02:1,'.k');
ylabel('Classification rate')

ylim([0 1]); xlabel('LDA start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
box off

end