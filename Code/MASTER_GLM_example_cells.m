%% Which cell to load
example_cells = {'R057-2015-02-15-TT11-cell2.mat', 'R057-2015-02-26-TT11-cell1.mat',... % cue modality (inc, dec)
    'R053-2014-11-15-TT13-cell1.mat', 'R056-2015-06-05-TT7-cell1.mat',... % cue location (inc, dec)
     'R060-2015-01-03-TT7-cell1.mat', 'R053-2014-11-15-TT4-cell1.mat',... % cue outcome (inc, dec)
      'R053-2014-11-12-TT5-cell1.mat','R060-2015-01-04-TT15-cell1.mat'}; % multiple cue features (mod x outcome, loc x out, respectively) 'R053-2014-11-12-TT4-cell3.mat' - old modxout 'R053-2014-11-12-TT5-cell1.mat' - old modxloc
example_coding = [1 1 2 2 3 3 5 4]; % what to visualize, 1 = modality, 2 = location, 3 = outcome, 4 = loc x out, 5 = mod x outcome
example_fig_PETH = {9,10,11,12,25,26,27,28};
example_fig_raster = {[1 5],[2 6],[3 7],[4 8],[17 21],[18 22],[19 23],[20 24]};
example_types = {'Identity coding (increasing)','Identity coding (decreasing)',...
    'Location coding (increasing)','Location coding (decreasing)',...
    'Outcome coding (increasing)','Outcome coding (decreasing)',...
    'Identity & outcome coding','Location & outcome coding'};
figure

  for i = 1:length(example_cells)
    load(example_cells{i});
    meta = metadata;
disp(i);

switch example_coding(i)
    case 1
%% light v sound blocks
peak_value(1) = max(PETH.Trial.MEAN.trials_light_PETH(4001:7000));
peak_value(2) = max(PETH.Trial.MEAN.trials_sound_PETH(4001:7000));
maximum_value = max(peak_value);

    time_light = -5:.001:10;
    time_light = time_light(1:length(PETH.Trial.MEAN.trials_light_PETH));
    time_sound = -5:.001:10;
    time_sound = time_sound(1:length(PETH.Trial.MEAN.trials_sound_PETH));
    subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(time_light,PETH.Trial.MEAN.trials_light_PETH,PETH.Trial.SEM.trials_light_PETH,'-r',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(time_sound,PETH.Trial.MEAN.trials_sound_PETH,PETH.Trial.SEM.trials_sound_PETH,'-b',1);
    xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
    box off;
     xlabel('Time from cue onset (s)'); 
     y_values =[0 round(maximum_value/5)*5];
set(gca, 'Ytick',y_values);
     if i == 1
    y = ylabel('Firing rate (Hz)');
set(y, 'position', get(y,'position')+[.14,0,0]); 
     end
    set(gca,'FontSize',18)%,'YTick',[]);

    case 2
%% trials separated by arm location
peak_value(1) = max(PETH.Arm.MEAN.photosensor1(4001:7000));
peak_value(2) = max(PETH.Arm.MEAN.photosensor2(4001:7000));
peak_value(3) = max(PETH.Arm.MEAN.photosensor3(4001:7000));
peak_value(4) = max(PETH.Arm.MEAN.photosensor4(4001:7000));
maximum_value = max(peak_value);

    photosensor1_time = -5:.001:10;
    photosensor1_time = photosensor1_time(1:length(PETH.Arm.MEAN.photosensor1));
    photosensor2_time = -5:.001:10;
    photosensor2_time = photosensor2_time(1:length(PETH.Arm.MEAN.photosensor2));
    photosensor3_time = -5:.001:10;
    photosensor3_time = photosensor3_time(1:length(PETH.Arm.MEAN.photosensor3));
    photosensor4_time = -5:.001:10;
    photosensor4_time = photosensor4_time(1:length(PETH.Arm.MEAN.photosensor4));
    subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(photosensor1_time,PETH.Arm.MEAN.photosensor1,PETH.Arm.SEM.photosensor1,'-m',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(photosensor2_time,PETH.Arm.MEAN.photosensor2,PETH.Arm.SEM.photosensor2,'-b',1);
    shadedErrorBar(photosensor3_time,PETH.Arm.MEAN.photosensor3,PETH.Arm.SEM.photosensor3,'-k',1);
    shadedErrorBar(photosensor4_time,PETH.Arm.MEAN.photosensor4,PETH.Arm.SEM.photosensor4,'-y',1);
     xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
  box off;
xlabel('Time from cue onset (s)');
     y_values =[0 round(maximum_value/5)*5];
set(gca, 'Ytick',y_values);
set(gca,'FontSize',18)%,'YTick',[]);

    case 3    
%% rew v unrew trials
peak_value(1) = max(PETH.Trial.MEAN.trials_rew_PETH(4001:7000));
peak_value(2) = max(PETH.Trial.MEAN.trials_unrew_PETH(4001:7000));
maximum_value = max(peak_value);

    time_rew = -5:.001:10;
    time_rew = time_rew(1:length(PETH.Trial.MEAN.trials_rew_PETH));
    time_unrew = -5:.001:10;
    time_unrew = time_unrew(1:length(PETH.Trial.MEAN.trials_unrew_PETH));
    subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(time_rew,PETH.Trial.MEAN.trials_rew_PETH,PETH.Trial.SEM.trials_rew_PETH,'-r',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(time_unrew,PETH.Trial.MEAN.trials_unrew_PETH,PETH.Trial.SEM.trials_unrew_PETH,'-g',1);
    xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
    box off;
     xlabel('Time from cue onset (s)'); 
          y_values =[0 round(maximum_value/5)*5];
set(gca, 'Ytick',y_values);
     if i == 5
    y = ylabel('Firing rate (Hz)');
set(y, 'position', get(y,'position')+[.14,0,0]); 
     end
set(gca,'FontSize',18)%,'YTick',[]);

    case 4
%% loc x out
    peak_value(1) = max(PETH.Arm.MEAN.photosensor1_rew(4001:7000));
peak_value(2) = max(PETH.Arm.MEAN.photosensor2_rew(4001:7000));
peak_value(3) = max(PETH.Arm.MEAN.photosensor3_rew(4001:7000));
peak_value(4) = max(PETH.Arm.MEAN.photosensor4_rew(4001:7000));

    peak_value(5) = max(PETH.Arm.MEAN.photosensor1_unrew(4001:7000));
peak_value(6) = max(PETH.Arm.MEAN.photosensor2_unrew(4001:7000));
peak_value(7) = max(PETH.Arm.MEAN.photosensor3_unrew(4001:7000));
peak_value(8) = max(PETH.Arm.MEAN.photosensor4_unrew(4001:7000));
maximum_value = max(peak_value);

photosensor1_rew_time = -5:.001:10;
    photosensor1_rew_time = photosensor1_rew_time(1:length(PETH.Arm.MEAN.photosensor1_rew));
    photosensor2_rew_time = -5:.001:10;
    photosensor2_rew_time = photosensor2_rew_time(1:length(PETH.Arm.MEAN.photosensor2_rew));
    photosensor3_rew_time = -5:.001:10;
    photosensor3_rew_time = photosensor3_rew_time(1:length(PETH.Arm.MEAN.photosensor3_rew));
    photosensor4_rew_time = -5:.001:10;
    photosensor4_rew_time = photosensor4_rew_time(1:length(PETH.Arm.MEAN.photosensor4_rew));
    photosensor1_unrew_time = -5:.001:10;
    photosensor1_unrew_time = photosensor1_unrew_time(1:length(PETH.Arm.MEAN.photosensor1_unrew));
    photosensor2_unrew_time = -5:.001:10;
    photosensor2_unrew_time = photosensor2_unrew_time(1:length(PETH.Arm.MEAN.photosensor2_unrew));
    photosensor3_unrew_time = -5:.001:10;
    photosensor3_unrew_time = photosensor3_unrew_time(1:length(PETH.Arm.MEAN.photosensor3_unrew));
    photosensor4_unrew_time = -5:.001:10;
    photosensor4_unrew_time = photosensor4_unrew_time(1:length(PETH.Arm.MEAN.photosensor4_unrew));
    
    subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(photosensor1_rew_time,PETH.Arm.MEAN.photosensor1_rew,PETH.Arm.SEM.photosensor1_rew,'-m',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(photosensor2_rew_time,PETH.Arm.MEAN.photosensor2_rew,PETH.Arm.SEM.photosensor2_rew,'-b',1);
    shadedErrorBar(photosensor3_rew_time,PETH.Arm.MEAN.photosensor3_rew,PETH.Arm.SEM.photosensor3_rew,'-k',1);
    shadedErrorBar(photosensor4_rew_time,PETH.Arm.MEAN.photosensor4_rew,PETH.Arm.SEM.photosensor4_rew,'-y',1);
    shadedErrorBar(photosensor1_unrew_time,PETH.Arm.MEAN.photosensor1_unrew,PETH.Arm.SEM.photosensor1_unrew,'-r',1);
    shadedErrorBar(photosensor2_unrew_time,PETH.Arm.MEAN.photosensor2_unrew,PETH.Arm.SEM.photosensor2_unrew,'-c',1);
    shadedErrorBar(photosensor3_unrew_time,PETH.Arm.MEAN.photosensor3_unrew,PETH.Arm.SEM.photosensor3_unrew,'-k',1);
    shadedErrorBar(photosensor4_unrew_time,PETH.Arm.MEAN.photosensor4_unrew,PETH.Arm.SEM.photosensor4_unrew,'-g',1);
     xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
  box off;
xlabel('Time from cue onset (s)');
     y_values =[0 round(maximum_value/5)*5];
set(gca, 'Ytick',y_values);
set(gca,'FontSize',18)%,'YTick',[]);

    case 5
%% mod x out
peak_value(1) = max(PETH.Trial.MEAN.rew_trials_light(4001:7000));
peak_value(2) = max(PETH.Trial.MEAN.unrew_trials_light(4001:7000));
peak_value(3) = max(PETH.Trial.MEAN.rew_trials_sound(4001:7000));
peak_value(4) = max(PETH.Trial.MEAN.unrew_trials_sound(4001:7000));
maximum_value = max(peak_value);

    rew_time_light = -5:.001:10;
    rew_time_light = rew_time_light(1:length(PETH.Trial.MEAN.rew_trials_light));
    unrew_time_light = -5:.001:10;
    unrew_time_light = unrew_time_light(1:length(PETH.Trial.MEAN.unrew_trials_light));
    rew_time_sound = -5:.001:10;
    rew_time_sound = rew_time_sound(1:length(PETH.Trial.MEAN.rew_trials_sound));
    unrew_time_sound = -5:.001:10;
    unrew_time_sound = unrew_time_sound(1:length(PETH.Trial.MEAN.unrew_trials_sound));
    subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(rew_time_light,PETH.Trial.MEAN.rew_trials_light,PETH.Trial.SEM.rew_trials_light,'-r',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.4),'.','color','black');
    shadedErrorBar(unrew_time_light,PETH.Trial.MEAN.unrew_trials_light,PETH.Trial.SEM.unrew_trials_light,'-g',1);
    xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.4)]);
    shadedErrorBar(rew_time_sound,PETH.Trial.MEAN.rew_trials_sound,PETH.Trial.SEM.rew_trials_sound,'-b',1);
    shadedErrorBar(unrew_time_sound,PETH.Trial.MEAN.unrew_trials_sound,PETH.Trial.SEM.unrew_trials_sound,'-c',1);
      xlabel('Time from cue onset (s)'); 
           y_values =[0 round(maximum_value/5)*5];
set(gca, 'Ytick',y_values);
%     ylabel('Firing rate (Hz)');
 box off;
  set(gca,'FontSize',18)%,'YTick',[]);

end

switch strcmp(sesh.session_id(1:4),'R060')
    case 0
%% rasterplots
subtightplot(8,4,example_fig_raster{i},[0,.03]);
for j = 1:length(metadata.TrialInfo_block1.trialT)
%     switch metadata.TrialInfo_block1.approached(j)
%         case 1
%             colour = 'black';
%         case 0
%             colour = 'red';
%     end
switch sesh.block_order
    case 1
        switch metadata.TrialInfo_block1.rewarded(j)
            case 1
                colour = 'red';
            case 0
                colour = 'green';
        end
    case 2
        switch metadata.TrialInfo_block1.rewarded(j)
            case 1
                colour = 'blue';
            case 0
                colour = 'cyan';
        end
end    
%     line([0 metadata.TrialInfo_block1.trial_length_analysis(j)],[j-.5 j-.5],'LineWidth',2,'color','yellow');
    if metadata.TrialInfo_block1.nosepoke_length(j) > 0
        if metadata.TrialInfo_block1.nosepoke_length(j) < 1 %limit data analysis to first 6 secs from nosepoke
           line([metadata.TrialInfo_block1.trial_length_analysis(j) metadata.TrialInfo_block1.trial_length_analysis(j)+metadata.TrialInfo_block1.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','magenta'); 
        elseif metadata.TrialInfo_block1.nosepoke_length(j) < 2
           line([metadata.TrialInfo_block1.trial_length_analysis(j) metadata.TrialInfo_block1.trial_length_analysis(j)+1],[j-.5 j-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo_block1.trial_length_analysis(j)+1 metadata.TrialInfo_block1.trial_length_analysis(j)+metadata.TrialInfo_block1.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','cyan');
        else
           line([metadata.TrialInfo_block1.trial_length_analysis(j) metadata.TrialInfo_block1.trial_length_analysis(j)+1],[j-.5 j-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo_block1.trial_length_analysis(j)+1 metadata.TrialInfo_block1.trial_length_analysis(j)+2],[j-.5 j-.5],'LineWidth',2,'color','cyan');
           line([metadata.TrialInfo_block1.trial_length_analysis(j)+2 metadata.TrialInfo_block1.trial_length_analysis(j)+metadata.TrialInfo_block1.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','green');
        end
            else
        line([metadata.TrialInfo_block1.trial_length_analysis(j) metadata.TrialInfo_block1.trial_length_analysis(j)+2],[j-.5 j-.5],'LineWidth',2,'color','green'); 
    end    
    for ii = 1:length(RAST.Trial{j})
        if RAST.Trial{j}(ii) > metadata.TrialInfo_block1.summary(j,9)
            break
        end
        line([RAST.Trial{j}(ii) RAST.Trial{j}(ii)],[j-1 j],'color',colour)
    end
end

for j = 1:length(metadata.TrialInfo_block2.trialT)
%     switch metadata.TrialInfo_block2.approached(j)
%         case 1
%             colour = 'black';
%         case 0
%             colour = 'red';
%     end
    
switch sesh.block_order
    case 2
        switch metadata.TrialInfo_block2.rewarded(j)
            case 1
                colour = 'red';
            case 0
                colour = 'green';
        end
    case 1
        switch metadata.TrialInfo_block2.rewarded(j)
            case 1
                colour = 'blue';
            case 0
                colour = 'cyan';
        end
end    
%     line([0 metadata.TrialInfo_block2.trial_length_analysis(j)],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','yellow'); 
    if metadata.TrialInfo_block2.nosepoke_length(j) > 0
        if metadata.TrialInfo_block2.nosepoke_length(j) < 1 %limit data analysis to first 6 secs from nosepoke
           line([metadata.TrialInfo_block2.trial_length_analysis(j) metadata.TrialInfo_block2.trial_length_analysis(j)+metadata.TrialInfo_block2.nosepoke_length(j)],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','magenta'); 
        elseif metadata.TrialInfo_block2.nosepoke_length(j) < 2
           line([metadata.TrialInfo_block2.trial_length_analysis(j) metadata.TrialInfo_block2.trial_length_analysis(j)+1],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo_block2.trial_length_analysis(j)+1 metadata.TrialInfo_block2.trial_length_analysis(j)+metadata.TrialInfo_block2.nosepoke_length(j)],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','cyan');
        else
           line([metadata.TrialInfo_block2.trial_length_analysis(j) metadata.TrialInfo_block2.trial_length_analysis(j)+1],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo_block2.trial_length_analysis(j)+1 metadata.TrialInfo_block2.trial_length_analysis(j)+2],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','cyan');
           line([metadata.TrialInfo_block2.trial_length_analysis(j)+2 metadata.TrialInfo_block2.trial_length_analysis(j)+metadata.TrialInfo_block2.nosepoke_length(j)],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','green');
        end
        else
        line([metadata.TrialInfo_block2.trial_length_analysis(j) metadata.TrialInfo_block2.trial_length_analysis(j)+2],[j+length(metadata.TrialInfo_block1.trialT)+10-.5 j+length(metadata.TrialInfo_block1.trialT)+10-.5],'LineWidth',2,'color','green'); 
    end    
    for ii = 1:length(RAST.Trial{j + length(metadata.TrialInfo_block1.trialT)})
        if RAST.Trial{j + length(metadata.TrialInfo_block1.trialT)}(ii) > metadata.TrialInfo_block2.summary(j,9)
            break
        end
        line([RAST.Trial{j + length(metadata.TrialInfo_block1.trialT)}(ii) RAST.Trial{j + length(metadata.TrialInfo_block1.trialT)}(ii)],[j + length(metadata.TrialInfo_block1.trialT)+10-1 j + length(metadata.TrialInfo_block1.trialT)+10],'color',colour)
    end
end
hold on; plot(0,0:10:length(RAST.Trial)+10,'.k'); %plot(-3:.5:6,length(metadata.TrialInfo_block1.trialT),'.b');
 xlim([-1.15 2.15]); ylim([0 length(RAST.Trial)+10+1]); set(gca,'XTick',[],'YTick',[]);
 if i == 1 || i == 5
ylabel('Trial number')
y_values = [length(metadata.TrialInfo_block1.trialT)/2   length(metadata.TrialInfo_block1.trialT)+(length(metadata.TrialInfo_block2.trialT)/2)];
y_labels ={'Block 1' 'Block 2'};
set(gca, 'Ytick',y_values,'YTickLabel',y_labels,'YTickLabelRotation',90);
 end
set(gca,'FontSize',18);
title(example_types{i})

    case 1
%% rasterplots - for newer rats
subtightplot(8,4,example_fig_raster{i},[0,.03]);
for j = 1:length(metadata.TrialInfo{1,1}.trialT)
%     switch metadata.TrialInfo{1,1}.approached(j)
%         case 1
%             colour = 'black';
%         case 0
%             colour = 'red';
%     end
switch sesh.block_order
    case 1
        switch metadata.TrialInfo{1,1}.rewarded(j)
            case 1
                colour = 'red';
            case 0
                colour = 'green';
        end
    case 2
        switch metadata.TrialInfo{1,1}.rewarded(j)
            case 1
                colour = 'blue';
            case 0
                colour = 'cyan';
        end
end
%     line([0 metadata.TrialInfo{1,1}.trial_length_analysis(j)],[j-.5 j-.5],'LineWidth',2,'color','yellow');
    if metadata.TrialInfo{1,1}.nosepoke_length(j) > 0
        if metadata.TrialInfo{1,1}.nosepoke_length(j) < 1 %limit data analysis to first 6 secs from nosepoke
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j) metadata.TrialInfo{1,1}.trial_length_analysis(j)+metadata.TrialInfo{1,1}.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','magenta'); 
        elseif metadata.TrialInfo{1,1}.nosepoke_length(j) < 2
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j) metadata.TrialInfo{1,1}.trial_length_analysis(j)+1],[j-.5 j-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j)+1 metadata.TrialInfo{1,1}.trial_length_analysis(j)+metadata.TrialInfo{1,1}.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','cyan');
        else
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j) metadata.TrialInfo{1,1}.trial_length_analysis(j)+1],[j-.5 j-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j)+1 metadata.TrialInfo{1,1}.trial_length_analysis(j)+2],[j-.5 j-.5],'LineWidth',2,'color','cyan');
           line([metadata.TrialInfo{1,1}.trial_length_analysis(j)+2 metadata.TrialInfo{1,1}.trial_length_analysis(j)+metadata.TrialInfo{1,1}.nosepoke_length(j)],[j-.5 j-.5],'LineWidth',2,'color','green');
        end
    else
        line([metadata.TrialInfo{1,1}.trial_length_analysis(j) metadata.TrialInfo{1,1}.trial_length_analysis(j)+metadata.TrialInfo{1,1}.trial_length_analysis(j+1)],[j-.5 j-.5],'LineWidth',2,'color','green'); 
    end    
    for ii = 1:length(RAST.Trial{j})
        if RAST.Trial{j}(ii) > metadata.TrialInfo{1,1}.summary(j,9)
            break
        end
        line([RAST.Trial{j}(ii) RAST.Trial{j}(ii)],[j-1 j],'color',colour)
    end
end

for j = 1:length(metadata.TrialInfo{1,2}.trialT)
%     switch metadata.TrialInfo{1,2}.approached(j)
%         case 1
%             colour = 'black';
%         case 0
%             colour = 'red';
%     end
    
switch sesh.block_order
    case 2
        switch metadata.TrialInfo{1,2}.rewarded(j)
            case 1
                colour = 'red';
            case 0
                colour = 'green';
        end
    case 1
        switch metadata.TrialInfo{1,2}.rewarded(j)
            case 1
                colour = 'blue';
            case 0
                colour = 'cyan';
        end
end
%     line([0 metadata.TrialInfo{1,2}.trial_length_analysis(j)],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','yellow');
    if metadata.TrialInfo{1,2}.nosepoke_length(j) > 0
        if metadata.TrialInfo{1,2}.nosepoke_length(j) < 1 %limit data analysis to first 6 secs from nosepoke
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j) metadata.TrialInfo{1,2}.trial_length_analysis(j)+metadata.TrialInfo{1,2}.nosepoke_length(j)],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','magenta'); 
        elseif metadata.TrialInfo{1,2}.nosepoke_length(j) < 2
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j) metadata.TrialInfo{1,2}.trial_length_analysis(j)+1],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j)+1 metadata.TrialInfo{1,2}.trial_length_analysis(j)+metadata.TrialInfo{1,2}.nosepoke_length(j)],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','cyan');
        else
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j) metadata.TrialInfo{1,2}.trial_length_analysis(j)+1],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','magenta');
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j)+1 metadata.TrialInfo{1,2}.trial_length_analysis(j)+2],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','cyan');
           line([metadata.TrialInfo{1,2}.trial_length_analysis(j)+2 metadata.TrialInfo{1,2}.trial_length_analysis(j)+metadata.TrialInfo{1,2}.nosepoke_length(j)],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','green');
        end
        else
        line([metadata.TrialInfo{1,2}.trial_length_analysis(j) metadata.TrialInfo{1,2}.trial_length_analysis(j)+metadata.TrialInfo{1,2}.trial_length_analysis(j+1)],[j+length(metadata.TrialInfo{1,1}.trialT)+10-.5 j+length(metadata.TrialInfo{1,1}.trialT)+10-.5],'LineWidth',2,'color','green'); 
    end    
    for ii = 1:length(RAST.Trial{j + length(metadata.TrialInfo{1,1}.trialT)})
        if RAST.Trial{j + length(metadata.TrialInfo{1,1}.trialT)}(ii) > metadata.TrialInfo{1,2}.summary(j,9)
            break
        end
        line([RAST.Trial{j + length(metadata.TrialInfo{1,1}.trialT)}(ii) RAST.Trial{j + length(metadata.TrialInfo{1,1}.trialT)}(ii)],[j + length(metadata.TrialInfo{1,1}.trialT)+10-1 j + length(metadata.TrialInfo{1,1}.trialT)+10],'color',colour)
    end
end
hold on; plot(0,0:10:length(RAST.Trial)+10,'.k'); %plot(-3:.5:6,length(metadata.TrialInfo{1,1}.trialT),'.b');
xlim([-1.15 2.15]); ylim([0 length(RAST.Trial)+10+1]); set(gca,'XTick',[],'YTick',[]); 
 if i == 1 || i == 5
ylabel('Trial number')
y_values = [length(metadata.TrialInfo{1,1}.trialT)/2   length(metadata.TrialInfo{1,1}.trialT)+(length(metadata.TrialInfo{1,2}.trialT)/2)];
y_labels ={'Block 1' 'Block 2'};
set(gca, 'Ytick',y_values,'YTickLabel',y_labels,'YTickLabelRotation',90);
 end
set(gca,'FontSize',18);
title(example_types{i})
end

clearvars -except load_var mat_files i example_cells example_coding example_fig_PETH example_fig_raster example_types
i = i + 1;
  end
  
%   tic; print(gcf,'-depsc','E:\examples2.eps'); toc;