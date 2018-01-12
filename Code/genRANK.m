function RANK = genRANK(DIFF,post_firing_rates,pre_firing_rates)
%  function RANK = genRANK(DIFF,post_firing_rates,pre_firing_rates)
%
%  Takes an input of resampled data within a pre- and post-event time
%  window, separates it according to pre- and post-, calculates the average
%  difference across events for each set of shuffles, and sorts according
%  to size.
%
%  INPUTS:
%    cfg.numSHUFF = 1000; how many rounds of resampled data are present
%    for each event
%    SHUFF; shuffled data to be separated into pre- and post-event 
%    event_time; list of event times, used to separated into pre- and
%    post-event
%    interval_length; list of interval lengths, can either be a single length
%    or matched for each event in SHUFF
%
%  OUTPUTS:
%    DIFF; cell array with sorted average differences of resampled data for
%    a pre- and post-event time window
%%
if isnan(DIFF) == 1
    RANK = NaN;
else
for i = 1:length(post_firing_rates)
        frate_diff_trial(i) = post_firing_rates(i) - pre_firing_rates(i);
end

avg_frate_diff_trial = mean(frate_diff_trial);

if isempty(find(avg_frate_diff_trial > DIFF,1,'last'))
    RANK = 0;
else
    RANK = find(avg_frate_diff_trial > DIFF,1,'last');
end
end