load_vars = find(ALL_matrix(:,2) == 1);

mat_files = dir('*.mat');

%%

for kk = 12:13%1:length(load_vars)
    load(mat_files(load_vars(kk)).name);
    disp(kk)
    figure(kk);
    subplot(2,1,1)
    plot(ACG.xbin2,ACG.ac2);
    ylim([0 .2])
    subplot(2,1,2)
    
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
    photosensor4_time = photosensor4_time(1:length(PETH.Arm.MEAN.photosensor4));shadedErrorBar(photosensor1_time,PETH.Arm.MEAN.photosensor1,PETH.Arm.SEM.photosensor1,'-m',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(photosensor2_time,PETH.Arm.MEAN.photosensor2,PETH.Arm.SEM.photosensor2,'-b',1);
    shadedErrorBar(photosensor3_time,PETH.Arm.MEAN.photosensor3,PETH.Arm.SEM.photosensor3,'-k',1);
    shadedErrorBar(photosensor4_time,PETH.Arm.MEAN.photosensor4,PETH.Arm.SEM.photosensor4,'-y',1);
     xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
  box off;
xlabel('Time from cue onset (s)');
%      y_values =[0 round(maximum_value/5)*5];
% set(gca, 'Ytick',y_values);
set(gca,'FontSize',16)%,'YTick',[]);

saveas(gcf,cat(2,'E:\cell',num2str(kk),'.png'));
close
end

%%
load_vars = find(ALL_matrix_np(:,2) == 1);

mat_files = dir('*.mat');

%%

for kk = 1:length(load_vars)
    load(mat_files(load_vars(kk)).name);
    disp(kk)
    figure(kk);
    subplot(2,1,1)
    plot(ACG.xbin2,ACG.ac2);
    ylim([0 .2])
    subplot(2,1,2)
    
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
    photosensor4_time = photosensor4_time(1:length(PETH.Arm.MEAN.photosensor4));shadedErrorBar(photosensor1_time,PETH.Arm.MEAN.photosensor1,PETH.Arm.SEM.photosensor1,'-m',1);
    hold on; plot(0,0:maximum_value/20:maximum_value+(maximum_value*.2),'.','color','black');
    shadedErrorBar(photosensor2_time,PETH.Arm.MEAN.photosensor2,PETH.Arm.SEM.photosensor2,'-b',1);
    shadedErrorBar(photosensor3_time,PETH.Arm.MEAN.photosensor3,PETH.Arm.SEM.photosensor3,'-k',1);
    shadedErrorBar(photosensor4_time,PETH.Arm.MEAN.photosensor4,PETH.Arm.SEM.photosensor4,'-y',1);
     xlim([-1.15 2.15]); ylim([0 maximum_value+(maximum_value*.2)]);
  box off;
xlabel('Time from cue onset (s)');
%      y_values =[0 round(maximum_value/5)*5];
% set(gca, 'Ytick',y_values);
set(gca,'FontSize',16)%,'YTick',[]);

saveas(gcf,cat(2,'E:\NPcell',num2str(kk),'.png'));
close
end