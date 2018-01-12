function SHUFF = genSHUFF(cfg_in,S,iv_out)
%  function SHUFF = genSHUFF(cfg_in,S,iv_out)
%
%  Takes an input of timestamps, and shuffles the placement of events within
%  a specified interval
%
%  INPUTS:
%    cfg.numSHUFF = 1000; how many rounds of shuffled data are to be generated
%    cfg.method = 1; method of bootstrapping. Can be 1 or 2. '1' generates x amount of
%        events using rand within the interval. '2' resamples keeping ISIs intact, 
%        including interval start and end times
%    S; vector of timestamps
%    iv_out; list of intervals, generated with iv
% 
%  OUTPUTS:
%    SHUFF; cell array with specified rounds of bootstrapping for the number
%        of intervals

cfg_def.numSHUFF = 1000;
cfg_def.method = 1;

cfg = ProcessConfig(cfg_def,cfg_in);

SHUFF = [];

switch cfg.method % which method of bootstrapping to use
    case 1
        for i = 1:length(iv_out.tstart) % for the amount of intervals
            temp = S(S > iv_out.tstart(i) & S < iv_out.tend(i)); % find events within interval
            if isempty(temp) % if empty move to next interval
%                 disp(cat(2,'No events detected in interval#',num2str(i)))
                SHUFF{i} = [];
                continue
            end
                        
            SHUFF{i} = (iv_out.tend(i)-iv_out.tstart(i)).*rand(length(temp),cfg.numSHUFF) + iv_out.tstart(i); % uses rand to place "temp" amount of events in interval
        end
        
    case 2
        for i = 1:length(iv_out.tstart) % for the amount of intervals
            temp = S(S > iv_out.tstart(i) & S < iv_out.tend(i)); % find events within interval
            if isempty(temp) % if empty move to next interval
%                 disp(cat(2,'No events detected in interval#',num2str(i)))
                SHUFF{i} = [];
                continue
            end
            
%             if length(temp) == 1
%                 disp(cat(2,'Warning: Interval#',num2str(i),' contains a single event'))
%             end
            
            S_shuff = [];
            S_shuff(1) = iv_out.tstart(i);
            S_shuff(2:length(temp)+1) = temp;
            S_shuff(length(temp)+2) = iv_out.tend(i);
            
            ISI = diff(S_shuff); % find ISIs
            
            for kk = 1:cfg.numSHUFF % for number of rounds of resampling
                idx = randperm(length(ISI)); % shuffle ISIs
                shuff_ISI = ISI(idx);
                SHUFF{i}(:,kk) = iv_out.tstart(i) + cumsum(shuff_ISI(1:end-1)); % generate resampled data
            end
        end
end