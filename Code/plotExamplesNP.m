function examplePLOTSnp = plotExamplesNP(directory)
% function examplePLOTSnp = plotExamplesNP(directory)
%
%
% INPUTS:
%
% OUTPUTS:

cd(directory)

%% Which cell to load
example_cells = [229 179 438 86 116 194];
example_coding = [1 1 2 2 3 3]; % what to visualize, 1 = modality, 2 = location, 3 = outcome
example_fig_PETH = {9,25,10,26,11,27};
example_fig_raster = {[1 5],[17 21],[2 6],[18 22],[3 7],[19 23]};
example_types = {'Identity coding (same)','Identity coding (different)',...
    'Location coding (same)','Location coding (different)',...
    'Outcome coding (same)','Outcome coding (different)',...
    'Identity & outcome coding','Identity & location coding'};
figure('units','normalized','outerposition',[0 0 1 1]);
mat_files = dir('*.mat');
for i = 1:length(example_cells)
    load(mat_files(example_cells(i)).name);
    meta = metadata;
    disp(i);
    
    switch example_coding(i)
        case 1
            %% light v sound blocks
            peak_value(1) = max(PETH.Nosepoke.MEAN.nosepoke_light_PETH(4001:7000));
            peak_value(2) = max(PETH.Nosepoke.MEAN.nosepoke_sound_PETH(4001:7000));
            maximum_value = max(peak_value);
            
            time_light = -5:.001:10;
            time_light = time_light(1:length(PETH.Nosepoke.MEAN.nosepoke_light_PETH));
            time_sound = -5:.001:10;
            time_sound = time_sound(1:length(PETH.Nosepoke.MEAN.nosepoke_sound_PETH));
            subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(time_light,PETH.Nosepoke.MEAN.nosepoke_light_PETH,PETH.Nosepoke.SEM.nosepoke_light_PETH,'-r',1);
            hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black'); plot(1,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','red');
            shadedErrorBar(time_sound,PETH.Nosepoke.MEAN.nosepoke_sound_PETH,PETH.Nosepoke.SEM.nosepoke_sound_PETH,'-b',1);
            xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
            box off;
            xlabel('Time from nosepoke (s)');
            y_values =[0 round(maximum_value/5)*5];
            set(gca, 'Ytick',y_values);
            %      if i == 1
            y = ylabel('Firing rate (Hz)');
            set(y, 'position', get(y,'position')+[.34,0,0]);
            %      end
            set(gca,'FontSize',18)%,'YTick',[]);
            
        case 2
            %% trials separated by arm location
            peak_value(1) = max(PETH.Receptacle.MEAN.receptacle1(4001:7000));
            peak_value(2) = max(PETH.Receptacle.MEAN.receptacle2(4001:7000));
            peak_value(3) = max(PETH.Receptacle.MEAN.receptacle3(4001:7000));
            peak_value(4) = max(PETH.Receptacle.MEAN.receptacle4(4001:7000));
            maximum_value = max(peak_value);
            
            photosensor1_time = -5:.001:10;
            photosensor1_time = photosensor1_time(1:length(PETH.Receptacle.MEAN.receptacle1));
            photosensor2_time = -5:.001:10;
            photosensor2_time = photosensor2_time(1:length(PETH.Receptacle.MEAN.receptacle2));
            photosensor3_time = -5:.001:10;
            photosensor3_time = photosensor3_time(1:length(PETH.Receptacle.MEAN.receptacle3));
            photosensor4_time = -5:.001:10;
            photosensor4_time = photosensor4_time(1:length(PETH.Receptacle.MEAN.receptacle4));
            subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(photosensor1_time,PETH.Receptacle.MEAN.receptacle1,PETH.Receptacle.SEM.receptacle1,'-m',1);
            hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black'); plot(1,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','red');
            shadedErrorBar(photosensor2_time,PETH.Receptacle.MEAN.receptacle2,PETH.Receptacle.SEM.receptacle2,'-b',1);
            shadedErrorBar(photosensor3_time,PETH.Receptacle.MEAN.receptacle3,PETH.Receptacle.SEM.receptacle3,'-k',1);
            shadedErrorBar(photosensor4_time,PETH.Receptacle.MEAN.receptacle4,PETH.Receptacle.SEM.receptacle4,'-y',1);
            xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
            box off;
            xlabel('Time from nosepoke (s)');
            y_values =[0 round(maximum_value/5)*5];
            set(gca, 'Ytick',y_values);
            set(gca,'FontSize',18)%,'YTick',[]);
            
        case 3
            %% rew v unrew trials
            peak_value(1) = max(PETH.Nosepoke.MEAN.nosepoke_rew_PETH(4001:7000));
            peak_value(2) = max(PETH.Nosepoke.MEAN.nosepoke_unrew_PETH(4001:7000));
            maximum_value = max(peak_value);
            
            time_rew = -5:.001:10;
            time_rew = time_rew(1:length(PETH.Nosepoke.MEAN.nosepoke_rew_PETH));
            time_unrew = -5:.001:10;
            time_unrew = time_unrew(1:length(PETH.Nosepoke.MEAN.nosepoke_unrew_PETH));
            subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(time_rew,PETH.Nosepoke.MEAN.nosepoke_rew_PETH,PETH.Nosepoke.SEM.nosepoke_rew_PETH,'-r',1);
            hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black'); plot(1,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','red');
            shadedErrorBar(time_unrew,PETH.Nosepoke.MEAN.nosepoke_unrew_PETH,PETH.Nosepoke.SEM.nosepoke_unrew_PETH,'-g',1);
            xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
            box off;
            xlabel('Time from nosepoke (s)');
            y_values =[0 round(maximum_value/5)*5];
            set(gca, 'Ytick',y_values);
            set(gca,'FontSize',18)
            
        case 4
            %% mod x loc
            peak_value(1) = max(PETH.Receptacle.MEAN.receptacle1_light(4001:7000));
            peak_value(2) = max(PETH.Receptacle.MEAN.receptacle2_light(4001:7000));
            peak_value(3) = max(PETH.Receptacle.MEAN.receptacle3_light(4001:7000));
            peak_value(4) = max(PETH.Receptacle.MEAN.receptacle4_light(4001:7000));
            
            peak_value(5) = max(PETH.Receptacle.MEAN.receptacle1_sound(4001:7000));
            peak_value(6) = max(PETH.Receptacle.MEAN.receptacle2_sound(4001:7000));
            peak_value(7) = max(PETH.Receptacle.MEAN.receptacle3_sound(4001:7000));
            peak_value(8) = max(PETH.Receptacle.MEAN.receptacle4_sound(4001:7000));
            maximum_value = max(peak_value);
            
            photosensor1_light_time = -5:.001:10;
            photosensor1_light_time = photosensor1_light_time(1:length(PETH.Receptacle.MEAN.receptacle1_light));
            photosensor2_light_time = -5:.001:10;
            photosensor2_light_time = photosensor2_light_time(1:length(PETH.Receptacle.MEAN.receptacle2_light));
            photosensor3_light_time = -5:.001:10;
            photosensor3_light_time = photosensor3_light_time(1:length(PETH.Receptacle.MEAN.receptacle3_light));
            photosensor4_light_time = -5:.001:10;
            photosensor4_light_time = photosensor4_light_time(1:length(PETH.Receptacle.MEAN.receptacle4_light));
            photosensor1_sound_time = -5:.001:10;
            photosensor1_sound_time = photosensor1_sound_time(1:length(PETH.Receptacle.MEAN.receptacle1_sound));
            photosensor2_sound_time = -5:.001:10;
            photosensor2_sound_time = photosensor2_sound_time(1:length(PETH.Receptacle.MEAN.receptacle2_sound));
            photosensor3_sound_time = -5:.001:10;
            photosensor3_sound_time = photosensor3_sound_time(1:length(PETH.Receptacle.MEAN.receptacle3_sound));
            photosensor4_sound_time = -5:.001:10;
            photosensor4_sound_time = photosensor4_sound_time(1:length(PETH.Receptacle.MEAN.receptacle4_sound));
            
            subtightplot(8,4,example_fig_PETH{i},[0,.03]); shadedErrorBar(photosensor1_light_time,PETH.Receptacle.MEAN.receptacle1_light,PETH.Receptacle.SEM.receptacle1_light,'-m',1);
            hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black'); plot(1,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','red');
            shadedErrorBar(photosensor2_light_time,PETH.Receptacle.MEAN.receptacle2_light,PETH.Receptacle.SEM.receptacle2_light,'-b',1);
            shadedErrorBar(photosensor3_light_time,PETH.Receptacle.MEAN.receptacle3_light,PETH.Receptacle.SEM.receptacle3_light,'-k',1);
            shadedErrorBar(photosensor4_light_time,PETH.Receptacle.MEAN.receptacle4_light,PETH.Receptacle.SEM.receptacle4_light,'-y',1);
            shadedErrorBar(photosensor1_sound_time,PETH.Receptacle.MEAN.receptacle1_sound,PETH.Receptacle.SEM.receptacle1_sound,'-r',1);
            shadedErrorBar(photosensor2_sound_time,PETH.Receptacle.MEAN.receptacle2_sound,PETH.Receptacle.SEM.receptacle2_sound,'-c',1);
            shadedErrorBar(photosensor3_sound_time,PETH.Receptacle.MEAN.receptacle3_sound,PETH.Receptacle.SEM.receptacle3_sound,'-k',1);
            shadedErrorBar(photosensor4_sound_time,PETH.Receptacle.MEAN.receptacle4_sound,PETH.Receptacle.SEM.receptacle4_sound,'-g',1);
            xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
            box off;
            xlabel('Time from nosepoke (s)');
            y_values =[0 round(maximum_value/5)*5];
            set(gca, 'Ytick',y_values);
            set(gca,'FontSize',18);
            
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
            xlabel('Time from nosepoke (s)');
            y_values =[0 round(maximum_value/5)*5];
            set(gca, 'Ytick',y_values);
            box off;
            set(gca,'FontSize',18)
            
    end
    
    subtightplot(8,4,example_fig_raster{i},[0,.03]);
    nosepoke = 1;
    
    for j = 1:length(metadata.TrialInfo{1,1}.trialT)
        if meta.TrialInfo{1,1}.nosepoke_length(j) ~= 0
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
            
            for ii = 1:length(RAST.Nosepoke{nosepoke})
                if RAST.Nosepoke{nosepoke}(ii) > meta.TrialInfo{1,1}.nosepoke_length(j)
                    break
                end
                line([RAST.Nosepoke{nosepoke}(ii) RAST.Nosepoke{nosepoke}(ii)],[nosepoke-1 nosepoke],'color',colour)
            end
            nosepoke = nosepoke + 1;
        end
    end   
    
    for j = 1:length(metadata.TrialInfo{1,2}.trialT)
        if meta.TrialInfo{1,2}.nosepoke_length(j) ~= 0
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
            
            for ii = 1:length(RAST.Nosepoke{nosepoke})
                if RAST.Nosepoke{nosepoke}(ii) > meta.TrialInfo{1,2}.nosepoke_length(j)
                    break
                end
                line([RAST.Nosepoke{nosepoke}(ii) RAST.Nosepoke{nosepoke}(ii)],[nosepoke-1+7.5 nosepoke+7.5],'color',colour)
            end
            nosepoke = nosepoke + 1;
        end
    end
    hold on; plot(0,0:5:length(RAST.Nosepoke)+10,'.k'); plot(1,0:5:length(RAST.Nosepoke)+10,'.r'); %plot(-3:.5:6,length(metadata.TrialInfo{1,1}.trialT),'.b');
    xlim([-1.15 2.15]); ylim([0 length(RAST.Nosepoke)+10+1]); set(gca,'XTick',[],'YTick',[]);
    if i == 1  || i == 2
        ylabel('Trial number')
        y_values = [nosepoke*.3   nosepoke*.85];
        y_labels ={'Block 1' 'Block 2'};
        set(gca, 'Ytick',y_values,'YTickLabel',y_labels,'YTickLabelRotation',90);
    end
    set(gca,'FontSize',18);
    title(example_types{i})
    
    clearvars -except load_var mat_files i example_cells example_coding example_fig_PETH example_fig_raster example_types
    i = i + 1;
end

end