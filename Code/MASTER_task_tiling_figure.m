figure;
subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
% colorbar; 
caxis([-3 4]);
title('Light block');
ylabel('Unit number');
xlabel('Time from cue onset (s)');
set(gca,'FontSize',16);

subplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]); 
 title('Sound block');
set(gca,'FontSize',16);
 
 subplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.sound_2_light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
 title('Sound vs. Light block');
set(gca,'FontSize',16);

 subplot(3,4,4)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_mod.light_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
 title('Light control');
set(gca,'FontSize',16);
 

 %%
 subplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
% colorbar; 
caxis([-3 4]);
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
title('Arm 1 trials');
set(gca,'FontSize',16);

subplot(3,4,6)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 2 trials');
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
set(gca,'FontSize',16);

subplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.arm2_2_arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 2 vs. Arm 1 trials');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
 ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
set(gca,'FontSize',16);

subplot(3,4,8)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_loc.arm1_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 1 control');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]);
 ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
set(gca,'FontSize',16);

 %%
 subplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
% colorbar; 
title('Rewarded trials');
caxis([-3 4]);
xlabel('Time from cue onset (s)'); ylabel('Unit number');
set(gca,'FontSize',16);

subplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.unrew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Unrewarded trials');
% colorbar; 
caxis([-3 4]);
 % set(gca,'YTick',[]); 
 xlabel('Time from cue onset (s)'); ylabel('Unit number');
 set(gca,'FontSize',16);
 
subplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.unrew_2_rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Unrewarded vs. Rewarded trials');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]); 
 xlabel('Time from cue onset (s)'); ylabel('Unit number');
set(gca,'FontSize',16);

subplot(3,4,12)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH_out.rew_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Rewarded control');
% colorbar; 
caxis([-3 4]);
% set(gca,'YTick',[]); 
 xlabel('Time from cue onset (s)'); ylabel('Unit number');
set(gca,'FontSize',16);
 
%%
figure
subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.light.zscore); 
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
% colorbar; 
 caxis([-3 4]);
 ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
 % set(gca,'YTick',[]); 
  title('Light block');
  set(gca,'FontSize',16);
  
subplot(3,4,2)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
%  colorbar; 
  caxis([-3 4]);
 %  set(gca,'YTick',[]); 
   title('Sound block');
   set(gca,'FontSize',16);
   
subplot(3,4,3)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.sound_2_light.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]); 
  title('Sound vs. Light block');
  ylabel('Unit number');
xlabel('Time from cue onset (s)');
set(gca,'FontSize',16);

subplot(3,4,4)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_mod.light_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]); 
  title('Light control');
  ylabel('Unit number');
xlabel('Time from cue onset (s)');
set(gca,'FontSize',16);

%%
subplot(3,4,5)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm1.zscore); 
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 1 trials');
% colorbar; 
 caxis([-3 4]);
 % set(gca,'YTick',[]);
  ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
set(gca,'FontSize',16);

subplot(3,4,6)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 2 trials');
%  colorbar; 
  caxis([-3 4]);
  ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
%   set(gca,'YTick',[]);
set(gca,'FontSize',16);

subplot(3,4,7)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.arm2_2_arm1.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 2 vs. Arm 1 trials');
%  colorbar; 
 caxis([-3 4]);
 ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
%  set(gca,'YTick',[]);
set(gca,'FontSize',16);

subplot(3,4,8)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_loc.arm1_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Arm 1 control');
%  colorbar; 
 caxis([-3 4]);
 ylabel('Unit number');
xlabel('Time from cue onset (s)'); 
%  set(gca,'YTick',[]);
set(gca,'FontSize',16);

%%
subplot(3,4,9)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Rewarded trials');
% colorbar; 
 caxis([-3 4]);
 % set(gca,'YTick',[]); 
  xlabel('Time from cue onset (s)'); 
   ylabel('Unit number');
  set(gca,'FontSize',16);
  
subplot(3,4,10)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.unrew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
title('Unrewarded trials');
%  colorbar; 
  caxis([-3 4]);
%   set(gca,'YTick',[]); 
   xlabel('Time from cue onset (s)');
    ylabel('Unit number');
   set(gca,'FontSize',16);
   
subplot(3,4,11)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH.unrew_2_rew.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]);
  xlabel('Time from cue onset (s)');
   ylabel('Unit number');
  title('Unrewarded vs. Rewarded trials');
set(gca,'FontSize',16);  

subplot(3,4,12)
imagesc(-1:.001:3,1:length(dir('*.mat')),MINsortedPETH_out.rew_2nd_2_1st.zscore);
hold on; plot(0,0:10:length(dir('*.mat')),'.','color','black');
%  colorbar; 
 caxis([-3 4]);
%  set(gca,'YTick',[]);
  xlabel('Time from cue onset (s)');
   ylabel('Unit number');
  title('Rewarded control');
set(gca,'FontSize',16);  
  
 %%
 figure
 subplot(3,4,1)
imagesc(-1:.001:3,1:length(dir('*.mat')),sortedPETH.light.zscore);
colorbar; 
  caxis([-3 4]);