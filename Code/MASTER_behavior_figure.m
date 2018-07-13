%% first load relevant variables from MASTER_behavior and MASTER_behavior_example

%%
group = [1 2 3 4];
figure
% subplot(2,3,[2 5])
subplot(7,3,[2 5 8 ])
gscatter([1 1 1 1],BEHAV_summary.APP.MEAN.rew_trials_light,group,'r','xo+*')
hold on;
gscatter([2 2 2 2],BEHAV_summary.APP.MEAN.unrew_trials_light,group,'g','xo+*')
gscatter([3 3 3 3],BEHAV_summary.APP.MEAN.rew_trials_sound,group,'b','xo+*')
gscatter([4 4 4 4],BEHAV_summary.APP.MEAN.unrew_trials_sound,group,'c','xo+*')
plot([.85 1.15],[mean_rew_light mean_rew_light],'k')
plot([1.85 2.15],[mean_unrew_light mean_unrew_light],'k')
plot([2.85 3.15],[mean_rew_sound mean_rew_sound],'k')
plot([3.85 4.15],[mean_unrew_sound mean_unrew_sound],'k')

xlim([0 5]); ylim([0 1]); %title('Proportion Approached');
% set(gca,'XTickLength', [0 0]); 
set(gca,'XTickLabel',{'','L1+','L2-','S1+','S2-'});
ylabel('Proportion approached'); xlabel('Cue type');
box off;
h = gca;
h.XRuler.TickLength = 0; 
set(h,'FontSize',20);
    
% figure
%subplot(2,3,[3 6])
subplot(7,3,[3 6 9])
gscatter([1 1 1 1],BEHAV_summary.Length.MEAN.rew_trials_light,group,'r','xo+*')
hold on;
gscatter([2 2 2 2],BEHAV_summary.Length.MEAN.unrew_trials_light,group,'g','xo+*')
gscatter([3 3 3 3],BEHAV_summary.Length.MEAN.rew_trials_sound,group,'b','xo+*')
gscatter([4 4 4 4],BEHAV_summary.Length.MEAN.unrew_trials_sound,group,'c','xo+*')
plot([.85 1.15],[mean_rew_light2 mean_rew_light2],'k')
plot([1.85 2.15],[mean_unrew_light2 mean_unrew_light2],'k')
plot([2.85 3.15],[mean_rew_sound2 mean_rew_sound2],'k')
plot([3.85 4.15],[mean_unrew_sound2 mean_unrew_sound2],'k')

xlim([0 5]); ylim([0 3]); %title('Trial Length');
% set(gca,'TickLength', [0 0]); box off;
set(gca,'XTickLabel',{'','L1+','L2-','S1+','S2-'})
    ylabel('Trial length (s)');  xlabel('Cue type'); %'XTickLabelRotation',90,
box off;
h = gca;
h.XRuler.TickLength = 0;   
set(h,'FontSize',20);

% subplot(2,3,1)
subplot(7,3,[1 4])
plot(prop_light_app_1,'color','r'); title('Light block'); ylabel('Proportion approached');
ylim([0 1.1]); hold on; plot(prop_light_app_2,'color','g'); 
plot(42,0:.05:1,'.','color','black'); xlim([0 70]);
% set(gca,'XTick', []);
xlabel('Session number');
box off;

light_sig_first = find(p_light_1 < .05, 1,'first')
light_sig_first_last = find(p_light_1(light_sig_first + 1:end) > .05, 1,'first') + light_sig_first
light_sig_next = find(p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_first_last
light_sig_next_last = find(p_light_1(light_sig_next + 1:end) > .05, 1,'first') + light_sig_next
light_sig_final = find(p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_next_last
light_sig_final_last = length(p_light_1)
plot([light_sig_first light_sig_first_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([light_sig_next light_sig_next_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([light_sig_final light_sig_final_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
set(gca,'FontSize',20);

% subplot(2,3,4)
subplot(7,3,[10 13])
plot(prop_sound_app_1,'color','c'); title('Sound block'); ylim([0 1.1]); %xlabel('Session number'); 
hold on; plot(prop_sound_app_2,'color','b'); plot(42,0:.05:1,'.','color','black'); 
xlim([0 70]); xlabel('Session number'); ylabel('Proportion approached'); 
box off;

sound_sig_first = find(p_sound_1 < .05, 1,'first')
sound_sig_first_last = find(p_sound_1(sound_sig_first + 1:end) > .05, 1,'first') + sound_sig_first
sound_sig_next = find(p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_first_last
sound_sig_next_last = find(p_sound_1(sound_sig_next + 1:end) > .05, 1,'first') + sound_sig_next
sound_sig_final = find(p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_next_last
sound_sig_final_last = length(p_sound_1)
plot([sound_sig_first sound_sig_first_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([sound_sig_next sound_sig_next_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
plot([sound_sig_final sound_sig_final_last], [.05 .05], '-k', 'LineWidth',3,'color','r')
set(gca,'FontSize',20);