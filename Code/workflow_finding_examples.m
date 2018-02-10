mat_files = dir('*.mat');
cell_id = [];
id_count = 1;
%%
for iCell = 360:length(ALL_matrix)
           if ALL_matrix(iCell,2) == 1 && ALL_matrix(iCell,3) == 1
    load(mat_files(iCell).name);
  cell_id{id_count} = mat_files(iCell).name;  
  sesh.PETH.Trial = 0;
  sesh.PETH.Nosepoke = 0;
  sesh.PETH.Receptacle = 0;
  sesh.PETH.Trial_np = 0;
  PETHs{id_count} = genPETH(sesh,metadata,spk_t,dataPoint);
  disp(id_count)
  id_count = id_count + 1;
           end
end
%%
for jCell = 1:length(PETHs)

figure    

%% mod x loc
    peak_value(1) = max(PETHs{jCell}.Arm.MEAN.photosensor1_rew(4001:7000));
peak_value(2) = max(PETHs{jCell}.Arm.MEAN.photosensor2_rew(4001:7000));
peak_value(3) = max(PETHs{jCell}.Arm.MEAN.photosensor3_rew(4001:7000));
peak_value(4) = max(PETHs{jCell}.Arm.MEAN.photosensor4_rew(4001:7000));

    peak_value(5) = max(PETHs{jCell}.Arm.MEAN.photosensor1_unrew(4001:7000));
peak_value(6) = max(PETHs{jCell}.Arm.MEAN.photosensor2_unrew(4001:7000));
peak_value(7) = max(PETHs{jCell}.Arm.MEAN.photosensor3_unrew(4001:7000));
peak_value(8) = max(PETHs{jCell}.Arm.MEAN.photosensor4_unrew(4001:7000));
maximum_value = max(peak_value);

photosensor1_rew_time = -5:.001:10;
    photosensor1_rew_time = photosensor1_rew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor1_rew));
    photosensor2_rew_time = -5:.001:10;
    photosensor2_rew_time = photosensor2_rew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor2_rew));
    photosensor3_rew_time = -5:.001:10;
    photosensor3_rew_time = photosensor3_rew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor3_rew));
    photosensor4_rew_time = -5:.001:10;
    photosensor4_rew_time = photosensor4_rew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor4_rew));
    photosensor1_unrew_time = -5:.001:10;
    photosensor1_unrew_time = photosensor1_unrew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor1_unrew));
    photosensor2_unrew_time = -5:.001:10;
    photosensor2_unrew_time = photosensor2_unrew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor2_unrew));
    photosensor3_unrew_time = -5:.001:10;
    photosensor3_unrew_time = photosensor3_unrew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor3_unrew));
    photosensor4_unrew_time = -5:.001:10;
    photosensor4_unrew_time = photosensor4_unrew_time(1:length(PETHs{jCell}.Arm.MEAN.photosensor4_unrew));
    
%     subtightplot(8,4,1,[0,.03]); 
    shadedErrorBar(photosensor1_rew_time,PETHs{jCell}.Arm.MEAN.photosensor1_rew,PETHs{jCell}.Arm.SEM.photosensor1_rew,'-m',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(photosensor2_rew_time,PETHs{jCell}.Arm.MEAN.photosensor2_rew,PETHs{jCell}.Arm.SEM.photosensor2_rew,'-b',1);
    shadedErrorBar(photosensor3_rew_time,PETHs{jCell}.Arm.MEAN.photosensor3_rew,PETHs{jCell}.Arm.SEM.photosensor3_rew,'-k',1);
    shadedErrorBar(photosensor4_rew_time,PETHs{jCell}.Arm.MEAN.photosensor4_rew,PETHs{jCell}.Arm.SEM.photosensor4_rew,'-y',1);
    shadedErrorBar(photosensor1_unrew_time,PETHs{jCell}.Arm.MEAN.photosensor1_unrew,PETHs{jCell}.Arm.SEM.photosensor1_unrew,'-r',1);
    shadedErrorBar(photosensor2_unrew_time,PETHs{jCell}.Arm.MEAN.photosensor2_unrew,PETHs{jCell}.Arm.SEM.photosensor2_unrew,'-c',1);
    shadedErrorBar(photosensor3_unrew_time,PETHs{jCell}.Arm.MEAN.photosensor3_unrew,PETHs{jCell}.Arm.SEM.photosensor3_unrew,'-k',1);
    shadedErrorBar(photosensor4_unrew_time,PETHs{jCell}.Arm.MEAN.photosensor4_unrew,PETHs{jCell}.Arm.SEM.photosensor4_unrew,'-g',1);
     xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
  box off;
xlabel('Time from cue onset (s)');
%      y_values =[0 round(maximum_value/5)*5];
% set(gca, 'Ytick',y_values);
set(gca,'FontSize',16)%,'YTick',[]);
           end
%            end

           %%