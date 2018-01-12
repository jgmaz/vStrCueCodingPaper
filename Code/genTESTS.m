function TESTS = genTESTS(FRATE)
% function TESTS = genTESTS(FRATE)
%
% statistical tests
%
% INPUTS:
% FRATE: struct of firing rates
%
% OUTPUTS:
% TESTS: MWU and WSR tests where appropriate

%% Wilxocon signed rank and Mann-whitney u-test
% WSR for task-related neurons
TESTS.WSR.Overall.Trial = signrank(FRATE.Task.Trial_firing_rate,FRATE.Overall.firing_rate_total);
TESTS.WSR.Task.Trial_b4_vs_Trial = signrank(FRATE.Task.Trial_B4_firing_rate,FRATE.Task.Trial_firing_rate); % pre-cue vs cue-onset

%MWU for trial type comparisons
TESTS.MWU.Cue.Trial = ranksum(FRATE.Cue.Trial_firing_rate_block1,FRATE.Cue.Trial_firing_rate_block2); % light v tone for cue-onset
TESTS.MWU.Reward.Trial = ranksum(FRATE.Reward.Trial_firing_rate_reward,FRATE.Reward.Trial_firing_rate_unreward); % rew v unrew for trial
TESTS.MWU.Approach.Trial = ranksum(FRATE.Approach.Trial_firing_rate_skip,FRATE.Approach.Trial_firing_rate_app); % app v skip for cue-onset

%KW for location comparisons
TESTS.KW.Arm.Trial = kruskalwallis(FRATE.Arm.Trial_firing_rate,FRATE.Arm.Trial_firing_rate_groups); % arm 1 v 2 v 3 v 4 for trial
close; close;