function Recodeplots = plotRecode(GLM_coeff,which_plot)
% function Recodeplots = plotRecode(GLM_coeff,which_plot)
%
%
% INPUTS:
%
% OUTPUTS:

Epoch = {'cueon' 'NP' 'outcome'}; %1 = cue on, 2 = NP, 3 = outcome
Predictors = {'Modality' 'Location' 'Outcome'};

rColorMap = [linspace(253/255, 255/255, 45),linspace(255/255, 49/255, 211)]; %77 253
gColorMap = [linspace(224/255, 255/255, 45),linspace(255/255, 163/255, 211)]; %146 224 49,163,84
bColorMap = [linspace(239/255, 255/255, 45),linspace(255/255, 84/255, 211)]; %33 239
colorMap = [rColorMap; gColorMap; bColorMap]';

mincolor = -.2;
maxcolor = .95;

%%
if which_plot == 2
    
    graph_title = {'Identity coding across task epochs' 'Location coding across task epochs' 'Outcome coding across task epochs'};
    
    colors.Modality = {[49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};
    colors.Location = {[49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};
    colors.Outcome = {[49/255 163/255 84/255] [255/255 207/255 250/255] [.8 .8 .8] ...
        [255/255 207/255 250/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [.8 .8 .8] [49/255 163/255 84/255] [49/255 163/255 84/255]};
    
    for iPred = 1:3
        figure
        heatmap(GLM_coeff.summary.(Predictors{iPred}).Corr,[],[],[],'ColorMap',colorMap,'Colorbar',true, ...
            'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [.6 .6 .6], 'TickAngle',90, 'ShowAllTicks', true)%...
        hold on
        
        plot([6.32 6.32],[0.65 6.3],'color',colors.(Predictors{iPred}){1},'LineWidth',10); plot([0.65 0.65],[0.65 6.3],'color',colors.(Predictors{iPred}){1},'LineWidth',10);
        plot([0.65 6.32],[6.32 6.3],'color',colors.(Predictors{iPred}){1},'LineWidth',10); plot([0.65 6.32],[0.65 0.65],'color',colors.(Predictors{iPred}){1},'LineWidth',10);
        plot([12.32 12.32],[0.65 6.3],'color',colors.(Predictors{iPred}){2},'LineWidth',10); plot([6.68 6.68],[0.65 6.3],'color',colors.(Predictors{iPred}){2},'LineWidth',10);
        plot([6.68 12.32],[6.3 6.3],'color',colors.(Predictors{iPred}){2},'LineWidth',10); plot([6.68 12.32],[0.65 0.65],'color',colors.(Predictors{iPred}){2},'LineWidth',10);
        plot([18.35 18.35],[0.65 6.3],'color',colors.(Predictors{iPred}){3},'LineWidth',10); plot([12.68 12.68],[0.65 6.3],'color',colors.(Predictors{iPred}){3},'LineWidth',10);
        plot([12.68 18.35],[6.3 6.3],'color',colors.(Predictors{iPred}){3},'LineWidth',10); plot([12.68 18.35],[0.65 0.65],'color',colors.(Predictors{iPred}){3},'LineWidth',10);
        
        plot([6.32 6.32],[6.7 12.3],'color',colors.(Predictors{iPred}){4},'LineWidth',10); plot([0.65 0.65],[6.7 12.3],'color',colors.(Predictors{iPred}){4},'LineWidth',10);
        plot([0.65 6.32],[12.3 12.3],'color',colors.(Predictors{iPred}){4},'LineWidth',10); plot([0.65 6.32],[6.7 6.7],'color',colors.(Predictors{iPred}){4},'LineWidth',10);
        plot([12.32 12.32],[6.7 12.3],'color',colors.(Predictors{iPred}){5},'LineWidth',10); plot([6.68 6.68],[6.7 12.3],'color',colors.(Predictors{iPred}){5},'LineWidth',10);
        plot([6.68 12.32],[12.3 12.3],'color',colors.(Predictors{iPred}){5},'LineWidth',10); plot([6.68 12.32],[6.7 6.7],'color',colors.(Predictors{iPred}){5},'LineWidth',10);
        plot([18.35 18.35],[6.7 12.3],'color',colors.(Predictors{iPred}){6},'LineWidth',10); plot([12.68 12.68],[6.7 12.3],'color',colors.(Predictors{iPred}){6},'LineWidth',10);
        plot([12.68 18.35],[12.3 12.3],'color',colors.(Predictors{iPred}){6},'LineWidth',10); plot([12.68 18.35],[6.7 6.7],'color',colors.(Predictors{iPred}){6},'LineWidth',10);
        
        plot([6.32 6.32],[12.7 18.35],'color',colors.(Predictors{iPred}){7},'LineWidth',10); plot([0.65 0.65],[12.7 18.35],'color',colors.(Predictors{iPred}){7},'LineWidth',10);
        plot([0.65 6.32],[18.35 18.35],'color',colors.(Predictors{iPred}){7},'LineWidth',10); plot([0.65 6.32],[12.7 12.7],'color',colors.(Predictors{iPred}){7},'LineWidth',10);
        plot([12.32 12.32],[12.7 18.35],'color',colors.(Predictors{iPred}){8},'LineWidth',10); plot([6.68 6.68],[12.7 18.35],'color',colors.(Predictors{iPred}){8},'LineWidth',10);
        plot([6.68 12.32],[18.35 18.35],'color',colors.(Predictors{iPred}){8},'LineWidth',10); plot([6.68 12.32],[12.7 12.7],'color',colors.(Predictors{iPred}){8},'LineWidth',10);
        plot([18.35 18.35],[12.7 18.35],'color',colors.(Predictors{iPred}){9},'LineWidth',10); plot([12.68 12.68],[12.7 18.35],'color',colors.(Predictors{iPred}){9},'LineWidth',10);
        plot([12.68 18.35],[18.35 18.35],'color',colors.(Predictors{iPred}){9},'LineWidth',10); plot([12.68 18.35],[12.7 12.7],'color',colors.(Predictors{iPred}){9},'LineWidth',10);
        
        plot([6.5 6.5],[0.56 19.49],'k','LineWidth',3); plot([-0.51 18.44],[6.5 6.5],'k','LineWidth',3);
        plot([12.5 12.5],[0.56 19.49],'k','LineWidth',3); plot([-0.51 18.44],[12.5 12.5],'k','LineWidth',3);
        xlim([0.5 18.5])
        ylim([0.5 18.5])
        set(gca,'FontSize',18)
        set(gcf,'Position', [10, 10, 1150, 950])
        y = ylabel('Outcome                   Nosepoke                   Cue-onset');
        x = xlabel('Cue-onset                        Nosepoke                         Outcome');
        ax = gca;
        ax.Clipping = 'off';
        title(graph_title{iPred})
    end
    
    %%
else
    switch which_plot
        case 1
            start_Epoch = 1;
            end_Epoch = 1;
        case 3
            start_Epoch = 2;
            end_Epoch = 3;
    end
    
    graph_title = {'Coding of cue features at cue-onset' 'Coding of cue features at nosepoke' 'Coding of cue features at outcome receipt'};
    
    colors.cueon = {[49/255 163/255 84/255] [49/255 163/255 84/255] [.8 .8 .8] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [.8 .8 .8] [49/255 163/255 84/255] [49/255 163/255 84/255]};
    colors.NP = {[49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255]};
    colors.outcome = {[49/255 163/255 84/255] [49/255 163/255 84/255] [49/255 163/255 84/255] ...
        [49/255 163/255 84/255] [49/255 163/255 84/255] [.8 .8 .8] ...
        [49/255 163/255 84/255] [.8 .8 .8] [49/255 163/255 84/255]};
    
    for iEpoch = start_Epoch:end_Epoch
        figure
        
        heatmap(GLM_coeff.summary.(Epoch{iEpoch}).Corr,[],[],[],'ColorMap',colorMap,'Colorbar',true, ...
            'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'NaNColor', [0.6 0.6 0.6], 'TickAngle', 45, 'ShowAllTicks', true)%...
        hold on
        
        plot([6.32 6.32],[0.65 6.3],'color',colors.(Epoch{iEpoch}){1},'LineWidth',10); plot([0.65 0.65],[0.65 6.3],'color',colors.(Epoch{iEpoch}){1},'LineWidth',10);
        plot([0.65 6.32],[6.32 6.3],'color',colors.(Epoch{iEpoch}){1},'LineWidth',10); plot([0.65 6.32],[0.65 0.65],'color',colors.(Epoch{iEpoch}){1},'LineWidth',10);
        plot([12.32 12.32],[0.65 6.3],'color',colors.(Epoch{iEpoch}){2},'LineWidth',10); plot([6.68 6.68],[0.65 6.3],'color',colors.(Epoch{iEpoch}){2},'LineWidth',10);
        plot([6.68 12.32],[6.3 6.3],'color',colors.(Epoch{iEpoch}){2},'LineWidth',10); plot([6.68 12.32],[0.65 0.65],'color',colors.(Epoch{iEpoch}){2},'LineWidth',10);
        plot([18.35 18.35],[0.65 6.3],'color',colors.(Epoch{iEpoch}){3},'LineWidth',10); plot([12.68 12.68],[0.65 6.3],'color',colors.(Epoch{iEpoch}){3},'LineWidth',10);
        plot([12.68 18.35],[6.3 6.3],'color',colors.(Epoch{iEpoch}){3},'LineWidth',10); plot([12.68 18.35],[0.65 0.65],'color',colors.(Epoch{iEpoch}){3},'LineWidth',10);
        
        plot([6.32 6.32],[6.7 12.3],'color',colors.(Epoch{iEpoch}){4},'LineWidth',10); plot([0.65 0.65],[6.7 12.3],'color',colors.(Epoch{iEpoch}){4},'LineWidth',10);
        plot([0.65 6.32],[12.3 12.3],'color',colors.(Epoch{iEpoch}){4},'LineWidth',10); plot([0.65 6.32],[6.7 6.7],'color',colors.(Epoch{iEpoch}){4},'LineWidth',10);
        plot([12.32 12.32],[6.7 12.3],'color',colors.(Epoch{iEpoch}){5},'LineWidth',10); plot([6.68 6.68],[6.7 12.3],'color',colors.(Epoch{iEpoch}){5},'LineWidth',10);
        plot([6.68 12.32],[12.3 12.3],'color',colors.(Epoch{iEpoch}){5},'LineWidth',10); plot([6.68 12.32],[6.7 6.7],'color',colors.(Epoch{iEpoch}){5},'LineWidth',10);
        plot([18.35 18.35],[6.7 12.3],'color',colors.(Epoch{iEpoch}){6},'LineWidth',10); plot([12.68 12.68],[6.7 12.3],'color',colors.(Epoch{iEpoch}){6},'LineWidth',10);
        plot([12.68 18.35],[12.3 12.3],'color',colors.(Epoch{iEpoch}){6},'LineWidth',10); plot([12.68 18.35],[6.7 6.7],'color',colors.(Epoch{iEpoch}){6},'LineWidth',10);
        
        plot([6.32 6.32],[12.7 18.35],'color',colors.(Epoch{iEpoch}){7},'LineWidth',10); plot([0.65 0.65],[12.7 18.35],'color',colors.(Epoch{iEpoch}){7},'LineWidth',10);
        plot([0.65 6.32],[18.35 18.35],'color',colors.(Epoch{iEpoch}){7},'LineWidth',10); plot([0.65 6.32],[12.7 12.7],'color',colors.(Epoch{iEpoch}){7},'LineWidth',10);
        plot([12.32 12.32],[12.7 18.35],'color',colors.(Epoch{iEpoch}){8},'LineWidth',10); plot([6.68 6.68],[12.7 18.35],'color',colors.(Epoch{iEpoch}){8},'LineWidth',10);
        plot([6.68 12.32],[18.35 18.35],'color',colors.(Epoch{iEpoch}){8},'LineWidth',10); plot([6.68 12.32],[12.7 12.7],'color',colors.(Epoch{iEpoch}){8},'LineWidth',10);
        plot([18.35 18.35],[12.7 18.35],'color',colors.(Epoch{iEpoch}){9},'LineWidth',10); plot([12.68 12.68],[12.7 18.35],'color',colors.(Epoch{iEpoch}){9},'LineWidth',10);
        plot([12.68 18.35],[18.35 18.35],'color',colors.(Epoch{iEpoch}){9},'LineWidth',10); plot([12.68 18.35],[12.7 12.7],'color',colors.(Epoch{iEpoch}){9},'LineWidth',10);
        
        plot([6.5 6.5],[0.56 19.49],'k','LineWidth',3); plot([-0.51 18.44],[6.5 6.5],'k','LineWidth',3);
        plot([12.5 12.5],[0.56 19.49],'k','LineWidth',3); plot([-0.51 18.44],[12.5 12.5],'k','LineWidth',3);
        
        xlim([0.5 18.5])
        ylim([0.5 18.5])
        set(gca,'FontSize',18)
        set(gcf,'Position', [10, 10, 1150, 950])
        y = ylabel('Outcome coding          Location coding           Identity coding');
        x = xlabel('Identity coding               Location coding              Outcome coding');
        ax = gca;
        ax.Clipping = 'off';
        title(graph_title{iEpoch})
    end
end

end