function DIFF = genDIFF(cfg_in,SHUFF,event_time,interval_length)
%  function DIFF = genDIFF(cfg.in,SHUFF,event_time,interval_length)
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
cfg_def.numSHUFF = 1000;

cfg = ProcessConfig(cfg_def,cfg_in);

count_pre = [];
count_post = [];
frate_pre = [];
frate_post = [];
frate_diff = [];
avg_frate_diff = [];
       
for i = 1:length(SHUFF) % for number of events
    if isempty(SHUFF{i}) % if no shuffled events detected, differences are automatically 0
        frate_diff(i,1) = 0;
    else
        if length(interval_length) == 1 % make sure 'interval_length' is a suitable length
            int_length = interval_length;
        elseif length(interval_length) == length(SHUFF)
            int_length = interval_length(i);
        else
            error('Error: interval_length must be equal in length to SHUFF (or a single interval)')
        end
        for kl = 1:length(SHUFF{i}) % for rounds of shuffled data presented
            count_pre{i}(:,kl) = sum(SHUFF{i}(:,kl) < event_time(i)); % count cells in pre-event window
            count_post{i}(:,kl) = sum(SHUFF{i}(:,kl) > event_time(i));
            frate_pre{i}(1,kl) =  count_pre{i}(1,kl) / int_length; % firing rate in pre-event window
            frate_post{i}(1,kl) =  count_post{i}(1,kl) / int_length;
            frate_diff(i,kl) = frate_post{i}(1,kl) - frate_pre{i}(1,kl); % firing rate difference between pre- and post-event windows
        end
    end
end

if length(frate_diff(1,:)) == 1 % if there are no shuffles essentially
    DIFF = NaN;
else
for kl = 1:cfg.numSHUFF % for each round of shuffled data
    avg_frate_diff(kl) = mean(frate_diff(:,kl)); % average firing rate difference
end

DIFF = sort(avg_frate_diff); % sorted average firing rate difference
end