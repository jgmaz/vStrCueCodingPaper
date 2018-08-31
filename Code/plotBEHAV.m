function BEHAVplot = plotBEHAV(directory)
% function BEHAVplot = plotBEHAV(directory)
%
%
% INPUTS:
%
% OUTPUTS:
%

load(strcat(directory,'Behavior_summary.mat'))

group = [1 2 3 4];
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(7,2,[2 4 6 ])
gscatter([1 1 1 1],BEHAV_summary.APP.MEAN.rew_trials_light,group,'r','xo+*',15)
hold on;
gscatter([2 2 2 2],BEHAV_summary.APP.MEAN.unrew_trials_light,group,'g','xo+*',15)
gscatter([3 3 3 3],BEHAV_summary.APP.MEAN.rew_trials_sound,group,'b','xo+*',15)
gscatter([4 4 4 4],BEHAV_summary.APP.MEAN.unrew_trials_sound,group,'c','xo+*',15)
plot([.85 1.15],[BEHAV_summary.mean_rew_light BEHAV_summary.mean_rew_light],'k')
plot([1.85 2.15],[BEHAV_summary.mean_unrew_light BEHAV_summary.mean_unrew_light],'k')
plot([2.85 3.15],[BEHAV_summary.mean_rew_sound BEHAV_summary.mean_rew_sound],'k')
plot([3.85 4.15],[BEHAV_summary.mean_unrew_sound BEHAV_summary.mean_unrew_sound],'k')

xlim([0 5]); ylim([0 1]);
set(gca,'XTickLabel',{'','L1+','L2-','S1+','S2-'});
ylabel('Proportion approached'); xlabel('Cue type');
box off;
h = gca;
set(h,'FontSize',20);
legend('off')

end