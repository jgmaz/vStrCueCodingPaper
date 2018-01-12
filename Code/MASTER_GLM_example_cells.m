%% Which cell to load
example_cells = {'R057-2015-02-15-TT11-cell2.mat', 'R057-2015-02-26-TT11-cell1.mat',... % cue modality
    'R053-2014-11-15-TT13-cell1.mat', 'R056-2015-06-05-TT7-cell1.mat',... % cue location
     'R060-2015-01-03-TT7-cell1.mat', 'R053-2014-11-15-TT4-cell1.mat',... % cue outcome
      'R053-2014-11-12-TT4-cell3.mat','R053-2014-11-12-TT5-cell1.mat', }; % multiple cue features (mod x loc, mod x outcome)
example_coding = [1 1 2 2 3 3 5 4]; % what to visualize, 1 = modality, 2 = location, 3 = outcome, 4 = mod x loc, 5 = mod x outcome
example_fig_PETH = {9,10,11,12,25,26,27,28};
example_fig_raster = {[1 5],[2 6],[3 7],[4 8],[17 21],[18 22],[19 23],[20 24]};
example_types = {'Cue modality - increasing','Cue modality - decreasing',...
    'Cue location - increasing','Cue location - decreasing',...
    'Cue outcome - increasing','Cue outcome - decreasing',...
    'Cue modality x cue location','Cue modality x cue outcome'};
figure

  for i = 1:length(example_cells)
    load(example_cells{i});
    meta = metadata;
disp(i);
%%
peak_value(1) = max(PETH.Trial.MEAN.rew_trials_light(4001:7000));
peak_value(2) = max(PETH.Trial.MEAN.unrew_trials_light(4001:7000));
peak_value(3) = max(PETH.Trial.MEAN.rew_trials_sound(4001:7000));
peak_value(4) = max(PETH.Trial.MEAN.unrew_trials_sound(4001:7000));

peak_value(5) = max(PETH.Arm.MEAN.photosensor1(4001:7000));
peak_value(6) = max(PETH.Arm.MEAN.photosensor2(4001:7000));
peak_value(7) = max(PETH.Arm.MEAN.photosensor3(4001:7000));
peak_value(8) = max(PETH.Arm.MEAN.photosensor4(4001:7000));

if example_coding(i) == 4
    peak_value(9) = max(PETH.Arm.MEAN.photosensor1_light(4001:7000));
peak_value(10) = max(PETH.Arm.MEAN.photosensor2_light(4001:7000));
peak_value(11) = max(PETH.Arm.MEAN.photosensor3_light(4001:7000));
peak_value(12) = max(PETH.Arm.MEAN.photosensor4_light(4001:7000));

    peak_value(13) = max(PETH.Arm.MEAN.photosensor1_sound(4001:7000));
peak_value(14) = max(PETH.Arm.MEAN.photosensor2_sound(4001:7000));
peak_value(15) = max(PETH.Arm.MEAN.photosensor3_sound(4001:7000));
peak_value(16) = max(PETH.Arm.MEAN.photosensor4_sound(4001:7000));

end

maximum_value = max(peak_value);

%% waveforms (and Rat name and isolation variables)
% x_coor = [1:32; 38:69; 75:106; 112:143];
% subtightplot(8,4,1); plot(x_coor', sesh.wav.mWV,'color','k'); %title(cat(2,sesh.session_id,':   tt ',num2str(sesh.tt_number),'   cell ',num2str(sesh.cell_number))) % box off
% set(gca,'XTick',[],'YTick',[]);  box off; %ylabel('Microvolts'); % ,'YTick',[]);
% % text(1,max(sesh.wav.mWV(:))+((max(sesh.wav.mWV(:)))/4),cat(2,'L-ratio: ',num2str(sesh.ClustQual.CluSep.Lratio),'     Iso Dist: ',num2str(sesh.ClustQual.CluSep.IsolationDist)));

switch example_coding(i)
    case 1
%% light v sound blocks
    time_light = -5:.001:10;
    time_light = time_light(1:length(PETH.Trial.MEAN.trials_light_PETH));
    time_sound = -5:.001:10;
    time_sound = time_sound(1:length(PETH.Trial.MEAN.trials_sound_PETH));
    subtightplot(8,4,example_fig_PETH{i}); shadedErrorBar(time_light,PETH.Trial.MEAN.trials_light_PETH,PETH.Trial.SEM.trials_light_PETH,'-r',1);
    hold on; plot(0,0:maximum_value/25:maximum_value+.5,'.','color','black');
    shadedErrorBar(time_sound,PETH.Trial.MEAN.trials_sound_PETH,PETH.Trial.SEM.trials_sound_PETH,'-b',1);
    xlim([-1 2]); ylim([0 maximum_value+.5]);
    box off;
     xlabel('Time from cue onset (s)'); 
     if i == 1
    ylabel('Firing rate (Hz)');
     end
    set(gca,'FontSize',16,'YTick',[]);

    case 2
%% trials separated by arm location
    photosensor1_time = -5:.001:10;
    photosensor1_time = photosensor1_time(1:length(PETH.Arm.MEAN.photosensor1));
    photosensor2_time = -5:.001:10;
    photosensor2_time = photosensor2_time(1:length(PETH.Arm.MEAN.photosensor2));
    photosensor3_time = -5:.001:10;
    photosensor3_time = photosensor3_time(1:length(PETH.Arm.MEAN.photosensor3));
    photosensor4_time = -5:.001:10;
    photosensor4_time = photosensor4_time(1:length(PETH.Arm.MEAN.photosensor4));
    subtightplot(8,4,example_fig_PETH{i}); shadedErrorBar(photosensor1_time,PETH.Arm.MEAN.photosensor1,PETH.Arm.SEM.photosensor1,'-m',1);
    hold on; plot(0,0:maximum_value/25:maximum_value+.5,'.','color','black');
    shadedErrorBar(photosensor2_time,PETH.Arm.MEAN.photosensor2,PETH.Arm.SEM.photosensor2,'-b',1);
    shadedErrorBar(photosensor3_time,PETH.Arm.MEAN.photosensor3,PETH.Arm.SEM.photosensor3,'-k',1);
    shadedErrorBar(photosensor4_time,PETH.Arm.MEAN.photosensor4,PETH.Arm.SEM.photosensor4,'-y',1);
     xlim([-1 2]); ylim([0 maximum_value+.5]);
  box off;
xlabel('Time from cue onset (s)');
set(gca,'FontSize',16,'YTick',[]);

    case 3    
%% rew v unrew trials
    time_rew = -5:.001:10;
    time_rew = time_rew(1:length(PETH.Trial.MEAN.trials_rew_PETH));
    time_unrew = -5:.001:10;
    time_unrew = time_unrew(1:length(PETH.Trial.MEAN.trials_unrew_PETH));
    subtightplot(8,4,example_fig_PETH{i}); shadedErrorBar(time_rew,PETH.Trial.MEAN.trials_rew_PETH,PETH.Trial.SEM.trials_rew_PETH,'-r',1);
    hold on; plot(0,0:maximum_value/25:maximum_value+.5,'.','color','black');
    shadedErrorBar(time_unrew,PETH.Trial.MEAN.trials_unrew_PETH,PETH.Trial.SEM.trials_unrew_PETH,'-g',1);
    xlim([-1 2]); ylim([0 maximum_value+.5]);
    box off;
     xlabel('Time from cue onset (s)'); 
     if i == 5
    ylabel('Firing rate (Hz)');
     end
set(gca,'FontSize',16,'YTick',[]);

    case 4
%% mod x loc
photosensor1_light_time = -5:.001:10;
    photosensor1_light_time = photosensor1_light_time(1:length(PETH.Arm.MEAN.photosensor1_light));
    photosensor2_light_time = -5:.001:10;
    photosensor2_light_time = photosensor2_light_time(1:length(PETH.Arm.MEAN.photosensor2_light));
    photosensor3_light_time = -5:.001:10;
    photosensor3_light_time = photosensor3_light_time(1:length(PETH.Arm.MEAN.photosensor3_light));
    photosensor4_light_time = -5:.001:10;
    photosensor4_light_time = photosensor4_light_time(1:length(PETH.Arm.MEAN.photosensor4_light));
    photosensor1_sound_time = -5:.001:10;
    photosensor1_sound_time = photosensor1_sound_time(1:length(PETH.Arm.MEAN.photosensor1_sound));
    photosensor2_sound_time = -5:.001:10;
    photosensor2_sound_time = photosensor2_sound_time(1:length(PETH.Arm.MEAN.photosensor2_sound));
    photosensor3_sound_time = -5:.001:10;
    photosensor3_sound_time = photosensor3_sound_time(1:length(PETH.Arm.MEAN.photosensor3_sound));
    photosensor4_sound_time = -5:.001:10;
    photosensor4_sound_time = photosensor4_sound_time(1:length(PETH.Arm.MEAN.photosensor4_sound));
    
    subtightplot(8,4,example_fig_PETH{i}); shadedErrorBar(photosensor1_light_time,PETH.Arm.MEAN.photosensor1_light,PETH.Arm.SEM.photosensor1_light,'-m',1);
    hold on; plot(0,0:maximum_value/25:maximum_value+.5,'.','color','black');
    shadedErrorBar(photosensor2_light_time,PETH.Arm.MEAN.photosensor2_light,PETH.Arm.SEM.photosensor2_light,'-b',1);
    shadedErrorBar(photosensor3_light_time,PETH.Arm.MEAN.photosensor3_light,PETH.Arm.SEM.photosensor3_light,'-k',1);
    shadedErrorBar(photosensor4_light_time,PETH.Arm.MEAN.photosensor4_light,PETH.Arm.SEM.photosensor4_light,'-y',1);
    shadedErrorBar(photosensor1_sound_time,PETH.Arm.MEAN.photosensor1_sound,PETH.Arm.SEM.photosensor1_sound,'-r',1);
    shadedErrorBar(photosensor2_sound_time,PETH.Arm.MEAN.photosensor2_sound,PETH.Arm.SEM.photosensor2_sound,'-c',1);
    shadedErrorBar(photosensor3_sound_time,PETH.Arm.MEAN.photosensor3_sound,PETH.Arm.SEM.photosensor3_sound,'-k',1);
    shadedErrorBar(photosensor4_sound_time,PETH.Arm.MEAN.photosensor4_sound,PETH.Arm.SEM.photosensor4_sound,'-g',1);
     xlim([-1 2]); ylim([0 maximum_value+.5]);
  box off;
xlabel('Time from cue onset (s)');
set(gca,'FontSize',16,'YTick',[]);

    case 5
%% mod x out
    rew_time_light = -5:.001:10;
    rew_time_light = rew_time_light(1:length(PETH.Trial.MEAN.rew_trials_light));
    unrew_time_light = -5:.001:10;
    unrew_time_light = unrew_time_light(1:length(PETH.Trial.MEAN.unrew_trials_light));
    rew_time_sound = -5:.001:10;
    rew_time_sound = rew_time_sound(1:length(PETH.Trial.MEAN.rew_trials_sound));
    unrew_time_sound = -5:.001:10;
    unrew_time_sound = unrew_time_sound(1:length(PETH.Trial.MEAN.unrew_trials_sound));
    subtightplot(8,4,example_fig_PETH{i}); shadedErrorBar(rew_time_light,PETH.Trial.MEAN.rew_trials_light,PETH.Trial.SEM.rew_trials_light,'-r',1);
    hold on; plot(0,0:maximum_value/25:maximum_value+.5,'.','color','black');
    shadedErrorBar(unrew_time_light,PETH.Trial.MEAN.unrew_trials_light,PETH.Trial.SEM.unrew_trials_light,'-g',1);
    xlim([-1 2]); ylim([0 maximum_value+.5]);
    shadedErrorBar(rew_time_sound,PETH.Trial.MEAN.rew_trials_sound,PETH.Trial.SEM.rew_trials_sound,'-b',1);
    shadedErrorBar(unrew_time_sound,PETH.Trial.MEAN.unrew_trials_sound,PETH.Trial.SEM.unrew_trials_sound,'-c',1);
      xlabel('Time from cue onset (s)'); 
%     ylabel('Firing rate (Hz)');
 box off;
  set(gca,'FontSize',16,'YTick',[]);

end

switch strcmp(sesh.session_id(1:4),'R060')
    case 0
%% rasterplots
subtightplot(8,4,example_fig_raster{i});
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
 xlim([-1 2]); ylim([0 length(RAST.Trial)+10+1]); set(gca,'XTick',[],'YTick',[]);
 if i == 1 || i == 5
ylabel('Trial number')
 end
set(gca,'FontSize',16);
title(example_types{i})

    case 1
%% rasterplots - for newer rats
subtightplot(8,4,example_fig_raster{i});
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
xlim([-1 2]); ylim([0 length(RAST.Trial)+10+1]); set(gca,'XTick',[],'YTick',[]); 
 if i == 1 || i == 5
ylabel('Trial number')
 end
set(gca,'FontSize',16);
title(example_types{i})
end

clearvars -except load_var mat_files i example_cells example_coding example_fig_PETH example_fig_raster example_types
i = i + 1;
  end