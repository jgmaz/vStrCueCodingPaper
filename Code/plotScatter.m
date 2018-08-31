function Scatter = plotScatter(spike_directory,directory)
% function Scatter = plotScatter(spike_directory,directory)
%
%
% INPUTS:
%
% OUTPUTS:

load(strcat(directory,'GLM_cueon_DATA_0.mat'),'ALL_matrix')
cd(spike_directory)

mat_files = dir('*.mat');

unspecific_number = 1;
light_number = 1;
sound_number = 1;

for jj = 1:length(ALL_matrix)
    load(mat_files(jj).name);
    disp(jj);
    
     if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
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

end