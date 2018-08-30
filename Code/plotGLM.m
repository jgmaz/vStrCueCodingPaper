function GLMplots = plotGLM(GLM_window,which_plot)
% function GLMplots = plotGLM(GLM_window,which_plot)
%
%
% INPUTS:
%
% OUTPUTS:

figure('units','normalized','outerposition',[0 0 1 1]);
switch which_plot
    
    %% Cue-onset cue features
    case 1
        subplot(2,3,1)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Modality,GLM_window.SHUFF.cueon.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Location,GLM_window.SHUFF.cueon.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Outcome,GLM_window.SHUFF.cueon.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        
        for iTime = 1:length(GLM_window.Table.cueon.Outcome.Zscore_recode)
            if GLM_window.Table.cueon.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.cueon.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.cueon.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylabel('Proportion of units')
        ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        %% Cue-onset all task variables
    case 2
        subplot(2,3,1)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Modality,GLM_window.SHUFF.cueon.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Location,GLM_window.SHUFF.cueon.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Outcome,GLM_window.SHUFF.cueon.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        
        for iTime = 1:length(GLM_window.Table.cueon.Outcome.Zscore_recode)
            if GLM_window.Table.cueon.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.cueon.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.cueon.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylabel('Proportion of units')
        ylim([0 .5]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,2)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Approach,GLM_window.SHUFF.cueon.prop.ZSTD.Approach,'--b',1);
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Latency,GLM_window.SHUFF.cueon.prop.ZSTD.Latency,'--r',1);
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(5,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(6,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.cueon.Outcome.Zscore_recode)
            if GLM_window.Table.cueon.Approach.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.cueon.Latency.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylim([0 .5]); %title('Behavioral measures');
        xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,3)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Trial,GLM_window.SHUFF.cueon.prop.ZSTD.Trial,'--b',1);
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.prop.MEAN.Previous,GLM_window.SHUFF.cueon.prop.ZSTD.Previous,'--r',1);
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(7,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.cueon.prop(8,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.cueon.Outcome.Zscore_recode)
            if GLM_window.Table.cueon.Trial.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.cueon.Previous.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color','r')
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylim([0 .5]);
        xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,4)
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(1,:),GLM_window.DATA.cueon.Rsquared.SEM(1,:),'-b',1);
        hold on
        plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(2,:),GLM_window.DATA.cueon.Rsquared.SEM(2,:),'-r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(3,:),GLM_window.DATA.cueon.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Modality,GLM_window.SHUFF.cueon.Rsquared.SEM.Modality,'--b',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Location,GLM_window.SHUFF.cueon.Rsquared.SEM.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Outcome,GLM_window.SHUFF.cueon.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        ylim([1 8]); ylabel('Improvement to R-Squared (%)'); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,5)
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(4,:),GLM_window.DATA.cueon.Rsquared.SEM(4,:),'-b',1);
        hold on
        plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(5,:),GLM_window.DATA.cueon.Rsquared.SEM(5,:),'-r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Approach,GLM_window.SHUFF.cueon.Rsquared.SEM.Approach,'--b',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Latency,GLM_window.SHUFF.cueon.Rsquared.SEM.Latency,'--r',1);
        ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        subplot(2,3,6)
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(6,:),GLM_window.DATA.cueon.Rsquared.SEM(6,:),'-b',1);
        hold on
        plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.cueon.Rsquared.MEAN(7,:),GLM_window.DATA.cueon.Rsquared.SEM(7,:),'-r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Trial,GLM_window.SHUFF.cueon.Rsquared.SEM.Trial,'--b',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.cueon.Rsquared.MEAN.Previous,GLM_window.SHUFF.cueon.Rsquared.SEM.Previous,'--r',1);
        ylim([1 8]); xlabel('GLM start relative to cue-onset (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        %% Nosepoke & outcome proportion of units
    case 3
        subplot(2,3,2)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Modality,GLM_window.SHUFF.NP.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Location,GLM_window.SHUFF.NP.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Outcome,GLM_window.SHUFF.NP.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.NP.Outcome.Zscore_recode)
            if GLM_window.Table.NP.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.NP.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.NP.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylabel('Proportion of units')
        ylim([0 .5]); xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,3)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Modality,GLM_window.SHUFF.outcome.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Location,GLM_window.SHUFF.outcome.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Outcome,GLM_window.SHUFF.outcome.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.outcome.Outcome.Zscore_recode)
            if GLM_window.Table.outcome.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.outcome.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.outcome.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylim([0 .5]); xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        %% Nosepoke & outcome with Rsquared
    case 4
        subplot(2,3,2)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Modality,GLM_window.SHUFF.NP.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Location,GLM_window.SHUFF.NP.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.prop.MEAN.Outcome,GLM_window.SHUFF.NP.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.NP.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.NP.Outcome.Zscore_recode)
            if GLM_window.Table.NP.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.NP.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.NP.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylabel('Proportion of units')
        ylim([0 .5]); xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,3)
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Modality,GLM_window.SHUFF.outcome.prop.ZSTD.Modality,'--b',1);
        
        hold on
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Location,GLM_window.SHUFF.outcome.prop.ZSTD.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.prop.MEAN.Outcome,GLM_window.SHUFF.outcome.prop.ZSTD.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(2,:),'-ob','MarkerFaceColor','b','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(3,:),'-or','MarkerFaceColor','r','MarkerSize',4)
        plot(-.5:.1:.5,GLM_window.DATA.outcome.prop(4,:),'-o','color',[.3718 .7176 .3612],'MarkerFaceColor',[.3718 .7176 .3612],'MarkerSize',4)
        for iTime = 1:length(GLM_window.Table.outcome.Outcome.Zscore_recode)
            if GLM_window.Table.outcome.Modality.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.03 .03], '-k', 'LineWidth',2,'color','b')
            end
            if GLM_window.Table.outcome.Location.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.02 .02], '-k', 'LineWidth',2,'color','r')
            end
            if GLM_window.Table.outcome.Outcome.Zscore_recode(iTime) == 1
                plot([-.65+(iTime*.1) -.55+(iTime*.1)],[.01 .01], '-k', 'LineWidth',2,'color',[.3718 .7176 .3612])
            end
        end
        plot(-.05,0.01:.01:.5,'.k'); plot(-.45,0.01:.01:.5,'.k');
        ylim([0 .5]); xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
        subplot(2,3,5)
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.NP.Rsquared.MEAN(1,:),GLM_window.DATA.NP.Rsquared.SEM(1,:),'-b',1);
        hold on
        plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.NP.Rsquared.MEAN(2,:),GLM_window.DATA.NP.Rsquared.SEM(2,:),'-r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.NP.Rsquared.MEAN(3,:),GLM_window.DATA.NP.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.Rsquared.MEAN.Modality,GLM_window.SHUFF.NP.Rsquared.SEM.Modality,'--b',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.Rsquared.MEAN.Location,GLM_window.SHUFF.NP.Rsquared.SEM.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.NP.Rsquared.MEAN.Outcome,GLM_window.SHUFF.NP.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        ylim([1 10]);  xlabel('GLM start relative to nosepoke (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        ylabel('Improvement to R-Squared (%)');
        box off
        subplot(2,3,6)
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.outcome.Rsquared.MEAN(1,:),GLM_window.DATA.outcome.Rsquared.SEM(1,:),'-b',1);
        hold on
        plot(-.05,1:.2:10,'.k'); plot(-.45,1:.2:10,'.k');
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.outcome.Rsquared.MEAN(2,:),GLM_window.DATA.outcome.Rsquared.SEM(2,:),'-r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.DATA.outcome.Rsquared.MEAN(3,:),GLM_window.DATA.outcome.Rsquared.SEM(3,:),{'-','color',[.3718 .7176 .3612]},1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.Rsquared.MEAN.Modality,GLM_window.SHUFF.outcome.Rsquared.SEM.Modality,'--b',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.Rsquared.MEAN.Location,GLM_window.SHUFF.outcome.Rsquared.SEM.Location,'--r',1);
        shadedErrorBar(-.5:.1:.5,GLM_window.SHUFF.outcome.Rsquared.MEAN.Outcome,GLM_window.SHUFF.outcome.Rsquared.SEM.Outcome,{'--','color',[.3718 .7176 .3612]},1);
        ylim([1 10]);  xlabel('GLM start relative to outcome (s)'); set(gca,'FontSize',20); xlim([-.5 .5]);
        box off
        
end

end