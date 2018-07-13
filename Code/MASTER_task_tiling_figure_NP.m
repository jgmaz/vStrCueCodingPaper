Conditions = {'light' 'sound' 'sound_2_light' 'rew' 'unrew' 'unrew_2_rew' 'arm1' 'arm2' 'arm2_2_arm1'};
for iCond = 1:length(Conditions)
    Exclude = [];
        Exclude = find(isnan(sortedPETH.(Conditions{iCond}).zscore(:,1)));
      sortedPETH.(Conditions{iCond}).zscore(Exclude,:) = []; 
       Exclude = [];
        Exclude = find(isnan(MINsortedPETH.(Conditions{iCond}).zscore(:,1)));
      MINsortedPETH.(Conditions{iCond}).zscore(Exclude,:) = []; 
end

%%
figure;
subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% colorbar; 
caxis([-3 4]);
title('Light block');
ylabel('Unit number');
% xlabel('Time from nosepoke (s)');
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

subplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]); 
 title('Sound block');
set(gca,'FontSize',18,'YTick',[],'XTick',[]);
 
 subplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound_2_light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
 title('Sound vs. Light block');
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

%  subplot(3,4,4)
% imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_mod.light_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% % ylabel('Unit number');
% % xlabel('Time from nosepoke (s)'); 
% % colorbar; 
% caxis([-3 4]);
% % set(gca,'YTick',[]);
%  title('Light control');
% set(gca,'FontSize',18,'YTick',[],'XTick',[]);
 

 %%
 subplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% colorbar; 
caxis([-3 4]);
ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
title('Arm 1');
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

subplot(3,4,6)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Arm 2');
% ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

subplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2_2_arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Arm 2 vs. Arm 1');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
%  ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

% subplot(3,4,8)
% imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_loc.receptacle1_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% title('Arm 1 control');
% % colorbar; 
% caxis([-3 4]);
% % set(gca,'YTick',[]);
% %  ylabel('Unit number');
% % xlabel('Time from nosepoke (s)'); 
% set(gca,'FontSize',18,'YTick',[],'XTick',[]);

 %%
 subplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% colorbar; 
title('Reward-available');
caxis([-3 4]);
xlabel('Time from nosepoke (s)'); 
ylabel('Unit number');
set(gca,'FontSize',18,'YTick',[])%,'XTick',[]);

subplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.unrew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Reward-unavailable');
% colorbar; 
caxis([-3 4]);
 % set(gca,'YTick',[]); 
  xlabel('Time from nosepoke (s)'); %ylabel('Unit number');
 set(gca,'FontSize',18,'YTick',[])%,'XTick',[]);
 
subplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.unrew_2_rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Reward-unavailable vs. -available');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]); 
  xlabel('Time from nosepoke (s)'); %ylabel('Unit number');
set(gca,'FontSize',18,'YTick',[])%,'XTick',[]);

% subplot(3,4,12)
% imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_out.rew_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% title('Reward-available control');
% % colorbar; 
% caxis([-3 4]);
% % set(gca,'YTick',[]); 
% %  xlabel('Time from nosepoke (s)'); ylabel('Unit number');
% set(gca,'FontSize',18,'YTick',[],'XTick',[]);
 
%%
figure
subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.light.zscore); 
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% colorbar; 
 caxis([-3 4]);
 ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
 % set(gca,'YTick',[]); 
  title('Light block');
  set(gca,'FontSize',18,'YTick',[],'XTick',[]);
  
subplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
%  colorbar; 
  caxis([-3 4]);
 %  set(gca,'YTick',[]); 
   title('Sound block');
   set(gca,'FontSize',18,'YTick',[],'XTick',[]);
   
subplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound_2_light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]); 
  title('Sound vs. Light block');
%   ylabel('Unit number');
% xlabel('Time from nosepoke (s)');
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

% subplot(3,4,4)
% imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_mod.light_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% %  colorbar; 
%  caxis([-3 4]);
% %  set(gca,'YTick',[]); 
%   title('Light control');
% %   ylabel('Unit number');
% % xlabel('Time from nosepoke (s)');
% set(gca,'FontSize',18,'YTick',[],'XTick',[]);

%%
subplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm1.zscore); 
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Arm 1');
% colorbar; 
 caxis([-3 4]);
 % set(gca,'YTick',[]);
  ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

subplot(3,4,6)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Arm 2');
%  colorbar; 
  caxis([-3 4]);
%   ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
%   set(gca,'YTick',[]);
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

subplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2_2_arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Arm 2 vs. Arm 1');
%  colorbar; 
 caxis([-3 4]);
%  ylabel('Unit number');
% xlabel('Time from nosepoke (s)'); 
%  set(gca,'YTick',[]);
set(gca,'FontSize',18,'YTick',[],'XTick',[]);

% subplot(3,4,8)
% imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_loc.receptacle1_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% title('Arm 1 control');
% %  colorbar; 
%  caxis([-3 4]);
% %  ylabel('Unit number');
% % xlabel('Time from nosepoke (s)'); 
% %  set(gca,'YTick',[]);
% set(gca,'FontSize',18,'YTick',[],'XTick',[]);

%%
subplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Reward-available');
% colorbar; 
 caxis([-3 4]);
 % set(gca,'YTick',[]); 
  xlabel('Time from nosepoke (s)'); 
   ylabel('Unit number');
  set(gca,'FontSize',18,'YTick',[]);
  
subplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.unrew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
title('Reward-unavailable');
%  colorbar; 
  caxis([-3 4]);
%   set(gca,'YTick',[]); 
   xlabel('Time from nosepoke (s)');
%     ylabel('Unit number');
   set(gca,'FontSize',18,'YTick',[]);
   
subplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.unrew_2_rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]);
  xlabel('Time from nosepoke (s)');
%    ylabel('Unit number');
  title('Reward-unavailable vs. -available');
set(gca,'FontSize',18,'YTick',[]);  

% subplot(3,4,12)
% imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_out.rew_2nd_2_1st.zscore);
% hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black'); plot(1,0:10:length(dir('*.mat')),'.','color','red');
% %  colorbar; 
%  caxis([-3 4]);
% %  set(gca,'YTick',[]);
%   xlabel('Time from nosepoke (s)');
% %    ylabel('Unit number');
%   title('Reward-available control');
% set(gca,'FontSize',18,'YTick',[]);  
  
 %%
 figure
 subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.light.zscore);
colorbar; 
  caxis([-3 4]);