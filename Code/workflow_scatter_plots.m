%% fig 3b from Day(2010) - cue onset for cue type
% cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R060\Mat'));

mat_files = dir('*.mat');

unspecific_number = 1;
light_number = 1;
sound_number = 1;

for jj = 1:length(ALL_matrix)%(dir('*.mat'))
    load(mat_files(jj).name);
    disp(jj);
    if RANK.two.Trial > 975 || RANK.two.Trial < 26
        if TESTS.MWU.Cue.Trial < .01
 block_drift.block1_length(jj) = length(FRATE.Cue.Trial_firing_rate_block1);
                block_drift.block1_half(jj) = round(block_drift.block1_length(jj) / 2);
                block_drift.b1_1st_avg(jj) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(jj)));
                block_drift.b1_2nd_avg(jj) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(jj)+1:end));
                block_drift.MWU_b1(jj) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(jj)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(jj)+1:end));
                
                block_drift.block2_length(jj) = length(FRATE.Cue.Trial_firing_rate_block2);
                block_drift.block2_half(jj) = round(block_drift.block2_length(jj) / 2);
                block_drift.b2_1st_avg(jj) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(jj)));
                block_drift.b2_2nd_avg(jj) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(jj)+1:end));
                block_drift.MWU_b2(jj) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(jj)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(jj)+1:end));
                
                switch block_drift.MWU_b1(jj) < .01 || block_drift.MWU_b2(jj) < .01
                    case 0
    if ALL_matrix(jj,1) == 1

            switch sesh.block_order
                case 1
                    if mean(FRATE.Cue.Trial_firing_rate_block1) > mean(FRATE.Cue.Trial_firing_rate_block2)
                        frates.cue.trial.light.MEAN.light(light_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                        frates.cue.trial.light.SEM.light(light_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
                        frates.cue.trial.light.MEAN.sound(light_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                        frates.cue.trial.light.SEM.sound(light_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                        light_number = light_number + 1;
                    else
                        frates.cue.trial.sound.MEAN.light(sound_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                        frates.cue.trial.sound.SEM.light(sound_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
                        frates.cue.trial.sound.MEAN.sound(sound_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                        frates.cue.trial.sound.SEM.sound(sound_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                        sound_number = sound_number + 1;
                    end
                case 2
                    if mean(FRATE.Cue.Trial_firing_rate_block2) > mean(FRATE.Cue.Trial_firing_rate_block1)
                        frates.cue.trial.light.MEAN.light(light_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                        frates.cue.trial.light.SEM.light(light_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                        frates.cue.trial.light.MEAN.sound(light_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                        frates.cue.trial.light.SEM.sound(light_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
                        light_number = light_number + 1;
                    else
                        frates.cue.trial.sound.MEAN.light(sound_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                        frates.cue.trial.sound.SEM.light(sound_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                        frates.cue.trial.sound.MEAN.sound(sound_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                        frates.cue.trial.sound.SEM.sound(sound_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
                        sound_number = sound_number + 1;
                    end
            end
        else
            switch sesh.block_order
                case 1
                    frates.cue.trial.unspecific.MEAN.light(unspecific_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                    frates.cue.trial.unspecific.SEM.light(unspecific_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
                    frates.cue.trial.unspecific.MEAN.sound(unspecific_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                    frates.cue.trial.unspecific.SEM.sound(unspecific_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                case 2
                    frates.cue.trial.unspecific.MEAN.light(unspecific_number) = mean(FRATE.Cue.Trial_firing_rate_block2);
                    frates.cue.trial.unspecific.SEM.light(unspecific_number) = std(FRATE.Cue.Trial_firing_rate_block2)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block2)));
                    frates.cue.trial.unspecific.MEAN.sound(unspecific_number) = mean(FRATE.Cue.Trial_firing_rate_block1);
                    frates.cue.trial.unspecific.SEM.sound(unspecific_number) = std(FRATE.Cue.Trial_firing_rate_block1)/(sqrt(length(FRATE.Cue.Trial_firing_rate_block1)));
            end
            unspecific_number = unspecific_number + 1;
        end
    end
end
end
end

% %%
% figure;
% errorbar(frates.cue.trial.MEAN.light,frates.cue.trial.MEAN.sound,frates.cue.trial.SEM.sound,'.','MarkerSize',10);
% hold on;
% herrorbar(frates.cue.trial.MEAN.light,frates.cue.trial.MEAN.sound,frates.cue.trial.SEM.light,'.');
% xlim([0 30]); ylim([0 30]);
% plot([0 30],[0 30]);
% 
% %%
% figure;
% errorbarxy(frates.cue.trial.MEAN.light,frates.cue.trial.MEAN.sound,frates.cue.trial.SEM.light,frates.cue.trial.SEM.sound,'black');

%%
frates.cue.trial.unspecific.MEAN.light = cat(2,R057_frates.cue.trial.unspecific.MEAN.light,R056_frates.cue.trial.unspecific.MEAN.light,R053_frates.cue.trial.unspecific.MEAN.light,R060_frates.cue.trial.unspecific.MEAN.light);
frates.cue.trial.unspecific.SEM.light = cat(2,R057_frates.cue.trial.unspecific.SEM.light,R056_frates.cue.trial.unspecific.SEM.light,R053_frates.cue.trial.unspecific.SEM.light,R060_frates.cue.trial.unspecific.SEM.light);
frates.cue.trial.unspecific.MEAN.sound = cat(2,R057_frates.cue.trial.unspecific.MEAN.sound,R056_frates.cue.trial.unspecific.MEAN.sound,R053_frates.cue.trial.unspecific.MEAN.sound,R060_frates.cue.trial.unspecific.MEAN.sound);
frates.cue.trial.unspecific.SEM.sound = cat(2,R057_frates.cue.trial.unspecific.SEM.sound,R056_frates.cue.trial.unspecific.SEM.sound,R053_frates.cue.trial.unspecific.SEM.sound,R060_frates.cue.trial.unspecific.SEM.sound);

frates.cue.trial.light.MEAN.light = cat(2,R057_frates.cue.trial.light.MEAN.light,R056_frates.cue.trial.light.MEAN.light,R053_frates.cue.trial.light.MEAN.light,R060_frates.cue.trial.light.MEAN.light);
frates.cue.trial.light.SEM.light = cat(2,R057_frates.cue.trial.light.SEM.light,R056_frates.cue.trial.light.SEM.light,R053_frates.cue.trial.light.SEM.light,R060_frates.cue.trial.light.SEM.light);
frates.cue.trial.light.MEAN.sound = cat(2,R057_frates.cue.trial.light.MEAN.sound,R056_frates.cue.trial.light.MEAN.sound,R053_frates.cue.trial.light.MEAN.sound,R060_frates.cue.trial.light.MEAN.sound);
frates.cue.trial.light.SEM.sound = cat(2,R057_frates.cue.trial.light.SEM.sound,R056_frates.cue.trial.light.SEM.sound,R053_frates.cue.trial.light.SEM.sound,R060_frates.cue.trial.light.SEM.sound);

frates.cue.trial.sound.MEAN.light = cat(2,R057_frates.cue.trial.sound.MEAN.light,R056_frates.cue.trial.sound.MEAN.light,R053_frates.cue.trial.sound.MEAN.light,R060_frates.cue.trial.sound.MEAN.light);
frates.cue.trial.sound.SEM.light = cat(2,R057_frates.cue.trial.sound.SEM.light,R056_frates.cue.trial.sound.SEM.light,R053_frates.cue.trial.sound.SEM.light,R060_frates.cue.trial.sound.SEM.light);
frates.cue.trial.sound.MEAN.sound = cat(2,R057_frates.cue.trial.sound.MEAN.sound,R056_frates.cue.trial.sound.MEAN.sound,R053_frates.cue.trial.sound.MEAN.sound,R060_frates.cue.trial.sound.MEAN.sound);
frates.cue.trial.sound.SEM.sound = cat(2,R057_frates.cue.trial.sound.SEM.sound,R056_frates.cue.trial.sound.SEM.sound,R053_frates.cue.trial.sound.SEM.sound,R060_frates.cue.trial.sound.SEM.sound);
%%
figure;
errorbarxy(frates.cue.trial.unspecific.MEAN.light,frates.cue.trial.unspecific.MEAN.sound,frates.cue.trial.unspecific.SEM.light,frates.cue.trial.unspecific.SEM.sound,'black');
hold on;
errorbarxy(frates.cue.trial.light.MEAN.light,frates.cue.trial.light.MEAN.sound,frates.cue.trial.light.SEM.light,frates.cue.trial.light.SEM.sound,'r');
errorbarxy(frates.cue.trial.sound.MEAN.light,frates.cue.trial.sound.MEAN.sound,frates.cue.trial.sound.SEM.light,frates.cue.trial.sound.SEM.sound,'b');
xlim([0 30]); ylim([0 30]);
xlabel('Light firing rate (Hz)');
ylabel('Sound firing rate (Hz)');
title('Comparison across blocks');
plot([0 30],[0 30]);
set(gca,'FontSize',18);


%% fig 3b from Day(2010) - np for cue type
cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R057\Mat'));

mat_files = dir('*.mat');

unspecific_number = 1;
light_number = 1;
sound_number = 1;

for jj = 1:length(dir('*.mat'))
    load(mat_files(jj).name);
    
    if RANK.two.Nosepoke > 975 || RANK.two.Nosepoke < 26
        if TESTS.MWU.Cue.Nosepoke < .01
            switch sesh.block_order
                case 1
                    if mean(FRATE.Cue.Nosepoke_firing_rate_block1) > mean(FRATE.Cue.Nosepoke_firing_rate_block2)
                        frates.cue.np.light.MEAN.light(light_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                        frates.cue.np.light.SEM.light(light_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
                        frates.cue.np.light.MEAN.sound(light_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                        frates.cue.np.light.SEM.sound(light_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                        light_number = light_number + 1;
                    else
                        frates.cue.np.sound.MEAN.light(sound_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                        frates.cue.np.sound.SEM.light(sound_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
                        frates.cue.np.sound.MEAN.sound(sound_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                        frates.cue.np.sound.SEM.sound(sound_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                        sound_number = sound_number + 1;
                    end
                case 2
                    if mean(FRATE.Cue.Nosepoke_firing_rate_block2) > mean(FRATE.Cue.Nosepoke_firing_rate_block1)
                        frates.cue.np.light.MEAN.light(light_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                        frates.cue.np.light.SEM.light(light_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                        frates.cue.np.light.MEAN.sound(light_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                        frates.cue.np.light.SEM.sound(light_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
                        light_number = light_number + 1;
                    else
                        frates.cue.np.sound.MEAN.light(sound_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                        frates.cue.np.sound.SEM.light(sound_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                        frates.cue.np.sound.MEAN.sound(sound_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                        frates.cue.np.sound.SEM.sound(sound_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
                        sound_number = sound_number + 1;
                    end
            end
        else
            switch sesh.block_order
                case 1
                    frates.cue.np.unspecific.MEAN.light(unspecific_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                    frates.cue.np.unspecific.SEM.light(unspecific_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
                    frates.cue.np.unspecific.MEAN.sound(unspecific_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                    frates.cue.np.unspecific.SEM.sound(unspecific_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                case 2
                    frates.cue.np.unspecific.MEAN.light(unspecific_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block2);
                    frates.cue.np.unspecific.SEM.light(unspecific_number) = std(FRATE.Cue.Nosepoke_firing_rate_block2)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block2)));
                    frates.cue.np.unspecific.MEAN.sound(unspecific_number) = mean(FRATE.Cue.Nosepoke_firing_rate_block1);
                    frates.cue.np.unspecific.SEM.sound(unspecific_number) = std(FRATE.Cue.Nosepoke_firing_rate_block1)/(sqrt(length(FRATE.Cue.Nosepoke_firing_rate_block1)));
            end
            unspecific_number = unspecific_number + 1;
        end
    end
end

% %%
% figure;
% errorbar(frates.cue.np.MEAN.light,frates.cue.np.MEAN.sound,frates.cue.np.SEM.sound,'.','MarkerSize',10);
% hold on;
% herrorbar(frates.cue.np.MEAN.light,frates.cue.np.MEAN.sound,frates.cue.np.SEM.light,'.');
% xlim([0 30]); ylim([0 30]);
% plot([0 30],[0 30]);
% 
% %%
% figure;
% errorbarxy(frates.cue.np.MEAN.light,frates.cue.np.MEAN.sound,frates.cue.np.SEM.light,frates.cue.np.SEM.sound,'black');

%%
figure;
errorbarxy(frates.cue.np.unspecific.MEAN.light,frates.cue.np.unspecific.MEAN.sound,frates.cue.np.unspecific.SEM.light,frates.cue.np.unspecific.SEM.sound,'black');
hold on;
errorbarxy(frates.cue.np.light.MEAN.light,frates.cue.np.light.MEAN.sound,frates.cue.np.light.SEM.light,frates.cue.np.light.SEM.sound,'b');
errorbarxy(frates.cue.np.sound.MEAN.light,frates.cue.np.sound.MEAN.sound,frates.cue.np.sound.SEM.light,frates.cue.np.sound.SEM.sound,'g');
xlim([0 30]); ylim([0 30]);
xlabel('Light firing rate (Hz)');
ylabel('Tone firing rate (Hz)');
plot([0 30],[0 30]);


%% fig 3b from Day(2010) - cue onset for reward type
cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R060\Mat'));

mat_files = dir('*.mat');

unspecific_number = 1;
reward_number = 1;
unreward_number = 1;

for jj = 1:length(dir('*.mat'))
    load(mat_files(jj).name);
    disp(jj);
    if RANK.two.Trial > 975 || RANK.two.Trial < 26
        if TESTS.MWU.Reward.Trial < .01
            if mean(FRATE.Reward.Trial_firing_rate_reward) > mean(FRATE.Reward.Trial_firing_rate_unreward)
                frates.rew.trial.reward.MEAN.reward(reward_number) = mean(FRATE.Reward.Trial_firing_rate_reward);
                frates.rew.trial.reward.SEM.reward(reward_number) = std(FRATE.Reward.Trial_firing_rate_reward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_reward)));
                frates.rew.trial.reward.MEAN.unreward(reward_number) = mean(FRATE.Reward.Trial_firing_rate_unreward);
                frates.rew.trial.reward.SEM.unreward(reward_number) = std(FRATE.Reward.Trial_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_unreward)));
                reward_number = reward_number + 1;
            else
                frates.rew.trial.unreward.MEAN.reward(unreward_number) = mean(FRATE.Reward.Trial_firing_rate_reward);
                frates.rew.trial.unreward.SEM.reward(unreward_number) = std(FRATE.Reward.Trial_firing_rate_reward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_reward)));
                frates.rew.trial.unreward.MEAN.unreward(unreward_number) = mean(FRATE.Reward.Trial_firing_rate_unreward);
                frates.rew.trial.unreward.SEM.unreward(unreward_number) = std(FRATE.Reward.Trial_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_unreward)));
                unreward_number = unreward_number + 1;
            end
        else
            frates.rew.trial.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Trial_firing_rate_reward);
            frates.rew.trial.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Trial_firing_rate_reward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_reward)));
            frates.rew.trial.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Trial_firing_rate_unreward);
            frates.rew.trial.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Trial_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_unreward)));
            unspecific_number = unspecific_number + 1;
        end
%     else
%         frates.rew.trial.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Trial_firing_rate_reward);
%         frates.rew.trial.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Trial_firing_rate_reward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_reward)));
%         frates.rew.trial.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Trial_firing_rate_unreward);
%         frates.rew.trial.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Trial_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Trial_firing_rate_unreward)));
%         unspecific_number = unspecific_number + 1;
    end
end

%%
frates.rew.trial.unspecific.MEAN.reward = cat(2,R057_rew_frates.rew.trial.unspecific.MEAN.reward,R056_rew_frates.rew.trial.unspecific.MEAN.reward,R053_rew_frates.rew.trial.unspecific.MEAN.reward,R060_rew_frates.rew.trial.unspecific.MEAN.reward);
frates.rew.trial.unspecific.SEM.reward = cat(2,R057_rew_frates.rew.trial.unspecific.SEM.reward,R056_rew_frates.rew.trial.unspecific.SEM.reward,R053_rew_frates.rew.trial.unspecific.SEM.reward,R060_rew_frates.rew.trial.unspecific.SEM.reward);
frates.rew.trial.unspecific.MEAN.unreward = cat(2,R057_rew_frates.rew.trial.unspecific.MEAN.unreward,R056_rew_frates.rew.trial.unspecific.MEAN.unreward,R053_rew_frates.rew.trial.unspecific.MEAN.unreward,R060_rew_frates.rew.trial.unspecific.MEAN.unreward);
frates.rew.trial.unspecific.SEM.unreward = cat(2,R057_rew_frates.rew.trial.unspecific.SEM.unreward,R056_rew_frates.rew.trial.unspecific.SEM.unreward,R053_rew_frates.rew.trial.unspecific.SEM.unreward,R060_rew_frates.rew.trial.unspecific.SEM.unreward);

frates.rew.trial.reward.MEAN.reward = cat(2,R057_rew_frates.rew.trial.reward.MEAN.reward,R053_rew_frates.rew.trial.reward.MEAN.reward,R060_rew_frates.rew.trial.reward.MEAN.reward);
frates.rew.trial.reward.SEM.reward = cat(2,R057_rew_frates.rew.trial.reward.SEM.reward,R053_rew_frates.rew.trial.reward.SEM.reward,R060_rew_frates.rew.trial.reward.SEM.reward);
frates.rew.trial.reward.MEAN.unreward = cat(2,R057_rew_frates.rew.trial.reward.MEAN.unreward,R053_rew_frates.rew.trial.reward.MEAN.unreward,R060_rew_frates.rew.trial.reward.MEAN.unreward);
frates.rew.trial.reward.SEM.unreward = cat(2,R057_rew_frates.rew.trial.reward.SEM.unreward,R053_rew_frates.rew.trial.reward.SEM.unreward,R060_rew_frates.rew.trial.reward.SEM.unreward);

frates.rew.trial.unreward.MEAN.reward = cat(2,R057_rew_frates.rew.trial.unreward.MEAN.reward,R056_rew_frates.rew.trial.unreward.MEAN.reward,R053_rew_frates.rew.trial.unreward.MEAN.reward,R060_rew_frates.rew.trial.unreward.MEAN.reward);
frates.rew.trial.unreward.SEM.reward = cat(2,R057_rew_frates.rew.trial.unreward.SEM.reward,R056_rew_frates.rew.trial.unreward.SEM.reward,R053_rew_frates.rew.trial.unreward.SEM.reward,R060_rew_frates.rew.trial.unreward.SEM.reward);
frates.rew.trial.unreward.MEAN.unreward = cat(2,R057_rew_frates.rew.trial.unreward.MEAN.unreward,R056_rew_frates.rew.trial.unreward.MEAN.unreward,R053_rew_frates.rew.trial.unreward.MEAN.unreward,R060_rew_frates.rew.trial.unreward.MEAN.unreward);
frates.rew.trial.unreward.SEM.unreward = cat(2,R057_rew_frates.rew.trial.unreward.SEM.unreward,R056_rew_frates.rew.trial.unreward.SEM.unreward,R053_rew_frates.rew.trial.unreward.SEM.unreward,R060_rew_frates.rew.trial.unreward.SEM.unreward);


%%
figure;
errorbarxy(frates.rew.trial.unspecific.MEAN.reward,frates.rew.trial.unspecific.MEAN.unreward,frates.rew.trial.unspecific.SEM.reward,frates.rew.trial.unspecific.SEM.unreward,'black');
hold on;
errorbarxy(frates.rew.trial.reward.MEAN.reward,frates.rew.trial.reward.MEAN.unreward,frates.rew.trial.reward.SEM.reward,frates.rew.trial.reward.SEM.unreward,'r');
errorbarxy(frates.rew.trial.unreward.MEAN.reward,frates.rew.trial.unreward.MEAN.unreward,frates.rew.trial.unreward.SEM.reward,frates.rew.trial.unreward.SEM.unreward,'g');
xlim([0 15]); ylim([0 15]);
% xlabel('Reward firing rate (Hz)');
% ylabel('Unreward firing rate (Hz)');
plot([0 15],[0 15]);



%% fig 3b from Day(2010) - reward delivery for reward type
cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R057\Mat'));

mat_files = dir('*.mat');

unspecific_number = 1;
reward_number = 1;
unreward_number = 1;

for jj = 1:length(dir('*.mat'))
    load(mat_files(jj).name);
    
    if RANK.two.Click2cue > 975 || RANK.two.Click2cue < 26
        if TESTS.MWU.Reward.Click2cue < .01
            if mean(FRATE.Reward.Click2cue_firing_rate_reward) > mean(FRATE.Reward.Click2cue_firing_rate_unreward)
                frates.rew.rewdel.reward.MEAN.reward(reward_number) = mean(FRATE.Reward.Click2cue_firing_rate_reward);
                frates.rew.rewdel.reward.SEM.reward(reward_number) = std(FRATE.Reward.Click2cue_firing_rate_reward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_reward)));
                frates.rew.rewdel.reward.MEAN.unreward(reward_number) = mean(FRATE.Reward.Click2cue_firing_rate_unreward);
                frates.rew.rewdel.reward.SEM.unreward(reward_number) = std(FRATE.Reward.Click2cue_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_unreward)));
                reward_number = reward_number + 1;
            else
                frates.rew.rewdel.unreward.MEAN.reward(unreward_number) = mean(FRATE.Reward.Click2cue_firing_rate_reward);
                frates.rew.rewdel.unreward.SEM.reward(unreward_number) = std(FRATE.Reward.Click2cue_firing_rate_reward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_reward)));
                frates.rew.rewdel.unreward.MEAN.unreward(unreward_number) = mean(FRATE.Reward.Click2cue_firing_rate_unreward);
                frates.rew.rewdel.unreward.SEM.unreward(unreward_number) = std(FRATE.Reward.Click2cue_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_unreward)));
                unreward_number = unreward_number + 1;
            end
        else
            frates.rew.rewdel.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Click2cue_firing_rate_reward);
            frates.rew.rewdel.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Click2cue_firing_rate_reward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_reward)));
            frates.rew.rewdel.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Click2cue_firing_rate_unreward);
            frates.rew.rewdel.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Click2cue_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_unreward)));
            unspecific_number = unspecific_number + 1;
        end
%     else
%         frates.rew.rewdel.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Click2cue_firing_rate_reward);
%         frates.rew.rewdel.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Click2cue_firing_rate_reward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_reward)));
%         frates.rew.rewdel.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Click2cue_firing_rate_unreward);
%         frates.rew.rewdel.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Click2cue_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Click2cue_firing_rate_unreward)));
%         unspecific_number = unspecific_number + 1;
    end
end

%%
figure;
errorbarxy(frates.rew.rewdel.unspecific.MEAN.reward,frates.rew.rewdel.unspecific.MEAN.unreward,frates.rew.rewdel.unspecific.SEM.reward,frates.rew.rewdel.unspecific.SEM.unreward,'black');
hold on;
errorbarxy(frates.rew.rewdel.reward.MEAN.reward,frates.rew.rewdel.reward.MEAN.unreward,frates.rew.rewdel.reward.SEM.reward,frates.rew.rewdel.reward.SEM.unreward,'b');
errorbarxy(frates.rew.rewdel.unreward.MEAN.reward,frates.rew.rewdel.unreward.MEAN.unreward,frates.rew.rewdel.unreward.SEM.reward,frates.rew.rewdel.unreward.SEM.unreward,'g');
xlim([0 15]); ylim([0 15]);
xlabel('Reward firing rate (Hz)');
ylabel('Unreward firing rate (Hz)');
plot([0 15],[0 15]);


%% fig 3b from Day(2010) - np for reward type
cd(cat(2,'E:\Jimmie\Jimmie\Analysis\R057\Mat'));

mat_files = dir('*.mat');

unspecific_number = 1;
reward_number = 1;
unreward_number = 1;

for jj = 1:length(dir('*.mat'))
    load(mat_files(jj).name);
    
    if RANK.two.Nosepoke > 975 || RANK.two.Nosepoke < 26
        if TESTS.MWU.Reward.Nosepoke < .01
            if mean(FRATE.Reward.Nosepoke_firing_rate_reward) > mean(FRATE.Reward.Nosepoke_firing_rate_unreward)
                frates.rew.np.reward.MEAN.reward(reward_number) = mean(FRATE.Reward.Nosepoke_firing_rate_reward);
                frates.rew.np.reward.SEM.reward(reward_number) = std(FRATE.Reward.Nosepoke_firing_rate_reward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_reward)));
                frates.rew.np.reward.MEAN.unreward(reward_number) = mean(FRATE.Reward.Nosepoke_firing_rate_unreward);
                frates.rew.np.reward.SEM.unreward(reward_number) = std(FRATE.Reward.Nosepoke_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_unreward)));
                reward_number = reward_number + 1;
            else
                frates.rew.np.unreward.MEAN.reward(unreward_number) = mean(FRATE.Reward.Nosepoke_firing_rate_reward);
                frates.rew.np.unreward.SEM.reward(unreward_number) = std(FRATE.Reward.Nosepoke_firing_rate_reward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_reward)));
                frates.rew.np.unreward.MEAN.unreward(unreward_number) = mean(FRATE.Reward.Nosepoke_firing_rate_unreward);
                frates.rew.np.unreward.SEM.unreward(unreward_number) = std(FRATE.Reward.Nosepoke_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_unreward)));
                unreward_number = unreward_number + 1;
            end
        else
            frates.rew.np.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Nosepoke_firing_rate_reward);
            frates.rew.np.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Nosepoke_firing_rate_reward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_reward)));
            frates.rew.np.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Nosepoke_firing_rate_unreward);
            frates.rew.np.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Nosepoke_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_unreward)));
            unspecific_number = unspecific_number + 1;
        end
%     else
%         frates.rew.np.unspecific.MEAN.reward(unspecific_number) = mean(FRATE.Reward.Nosepoke_firing_rate_reward);
%         frates.rew.np.unspecific.SEM.reward(unspecific_number) = std(FRATE.Reward.Nosepoke_firing_rate_reward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_reward)));
%         frates.rew.np.unspecific.MEAN.unreward(unspecific_number) = mean(FRATE.Reward.Nosepoke_firing_rate_unreward);
%         frates.rew.np.unspecific.SEM.unreward(unspecific_number) = std(FRATE.Reward.Nosepoke_firing_rate_unreward)/(sqrt(length(FRATE.Reward.Nosepoke_firing_rate_unreward)));
%         unspecific_number = unspecific_number + 1;
    end
end

%%
figure;
errorbarxy(frates.rew.np.unspecific.MEAN.reward,frates.rew.np.unspecific.MEAN.unreward,frates.rew.np.unspecific.SEM.reward,frates.rew.np.unspecific.SEM.unreward,'black');
hold on;
errorbarxy(frates.rew.np.reward.MEAN.reward,frates.rew.np.reward.MEAN.unreward,frates.rew.np.reward.SEM.reward,frates.rew.np.reward.SEM.unreward,'b');
errorbarxy(frates.rew.np.unreward.MEAN.reward,frates.rew.np.unreward.MEAN.unreward,frates.rew.np.unreward.SEM.reward,frates.rew.np.unreward.SEM.unreward,'g');
xlim([0 15]); ylim([0 15]);
xlabel('Reward firing rate (Hz)');
ylabel('Unreward firing rate (Hz)');
plot([0 15],[0 15]);