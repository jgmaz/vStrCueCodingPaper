function BEHAVplot = plotBEHAV(behavior_directory,directory)
% function BEHAVplot = plotBEHAV(behavior_directory,directory)
%
%
% INPUTS:
%
% OUTPUTS:
%

load(strcat(directory,'Behavior_summary.mat'))
load(strcat(behavior_directory,'R060_learning_curve.mat'))

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

%%
subplot(7,2,[1 3])
plot(R060_learning_curve.prop_light_app_1,'color','r'); title('Light block'); ylabel('Proportion approached');
ylim([0 1.1]); hold on; plot(R060_learning_curve.prop_light_app_2,'color','g'); 
plot(42,0:.05:1,'.','color','black'); xlim([0 70]);
xlabel('Session number');
box off;

light_sig_first = find(R060_learning_curve.p_light_1 < .05, 1,'first');
light_sig_first_last = find(R060_learning_curve.p_light_1(light_sig_first + 1:end) > .05, 1,'first') + light_sig_first;
light_sig_next = find(R060_learning_curve.p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_first_last;
light_sig_next_last = find(R060_learning_curve.p_light_1(light_sig_next + 1:end) > .05, 1,'first') + light_sig_next;
light_sig_final = find(R060_learning_curve.p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_next_last;
light_sig_final_last = length(R060_learning_curve.p_light_1);
plot([light_sig_first light_sig_first_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([light_sig_next light_sig_next_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([light_sig_final light_sig_final_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
set(gca,'FontSize',20);

subplot(7,2,[7 9])
plot(R060_learning_curve.prop_sound_app_1,'color','c'); title('Sound block'); ylim([0 1.1]); %xlabel('Session number'); 
hold on; plot(R060_learning_curve.prop_sound_app_2,'color','b'); plot(42,0:.05:1,'.','color','black'); 
xlim([0 70]); xlabel('Session number'); ylabel('Proportion approached'); 
box off;

sound_sig_first = find(R060_learning_curve.p_sound_1 < .05, 1,'first');
sound_sig_first_last = find(R060_learning_curve.p_sound_1(sound_sig_first + 1:end) > .05, 1,'first') + sound_sig_first;
sound_sig_next = find(R060_learning_curve.p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_first_last;
sound_sig_next_last = find(R060_learning_curve.p_sound_1(sound_sig_next + 1:end) > .05, 1,'first') + sound_sig_next;
sound_sig_final = find(R060_learning_curve.p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_next_last;
sound_sig_final_last = length(R060_learning_curve.p_sound_1);
plot([sound_sig_first sound_sig_first_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([sound_sig_next sound_sig_next_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([sound_sig_final sound_sig_final_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
set(gca,'FontSize',20);
end