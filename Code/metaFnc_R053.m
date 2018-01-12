function meta = metaFnc(fname,block_order)
% function meta = metaFnc(fname,block_order)
%
% computes meta file combining behavioural and recording data
%
% INPUTS:
% fname: string
% block_order: 1 if light block came first, 2 if sound block came first
%
% OUTPUTS:
% meta.TrialInfo_block1: ...

disp('loading meta file')
switch block_order
    case 1
        behav_1 = '_light_behaviour.mat';
        behav_2 = '_sound_behaviour.mat';
        neural_1 = '-light.dat';
        neural_2 = '-sound.dat';
    case 2
        behav_1 = '_sound_behaviour.mat';
        behav_2 = '_light_behaviour.mat';
        neural_1 = '-sound.dat';
        neural_2 = '-light.dat';
end

% load behavioural data
load(cat(2,'C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\',fname(6:15),'\',fname,behav_1))
meta.TrialInfo_block1 = TrialFnc(eventLog); % generate behavioural data
eventLog_block1 = eventLog;
load(cat(2,'C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\',fname(6:15),'\',fname,behav_2))
meta.TrialInfo_block2 = TrialFnc(eventLog); % generate behavioural data
eventLog_block2 = eventLog;

% load AMPX maze data
fname_pre = cat(2,fname,'-pre.dat');
fname_block1 = cat(2,fname,neural_1);
fname_block2 = cat(2,fname,neural_2);
fname_post = cat(2,fname,'-post.dat');

channels_to_load = 65:76; %non neural signals
cd(cat(2,'E:\Jimmie\Data\',fname(1:4),'\',fname));
meta.data_pre = AMPX_loadData(fname_pre,channels_to_load);
meta.data_block1 = AMPX_loadData(fname_block1,channels_to_load);
meta.data_block2 = AMPX_loadData(fname_block2,channels_to_load);
meta.data_post = AMPX_loadData(fname_post,channels_to_load);

meta.gap1.data = zeros(((meta.data_block1.hdr.start_timestamp - meta.data_pre.hdr.end_timestamp)*20),1);
meta.gap2.data = zeros(((meta.data_block2.hdr.start_timestamp - meta.data_block1.hdr.end_timestamp)*20),1);
meta.gap3.data = zeros(((meta.data_post.hdr.start_timestamp - meta.data_block2.hdr.end_timestamp)*20),1);

meta.gap1.tvec(:,1) = (0:length(meta.gap1.data)-1) /20000;
meta.gap2.tvec(:,1) = (0:length(meta.gap2.data)-1) /20000;
meta.gap3.tvec(:,1) = (0:length(meta.gap3.data)-1) /20000;

%% Find time offset between matlab and AMPX for block 1
threshold = 3000; %threshold for an event
refractory_period = 5; %period of time before passing threshold considered a new event
trackPB_1 = find(meta.data_block1.channels{4}> threshold); %find possible events
trackPB_1_new = trackPB_1(1);
for i = 2:length(trackPB_1)
    if meta.data_block1.tvec(trackPB_1(i))- meta.data_block1.tvec(trackPB_1(i-1)) > refractory_period %new event if far away enough in time
        trackPB_1_new(end+1) = trackPB_1(i);
    end
end

trackPB_2 = find(meta.data_block1.channels{3}> threshold); %find possible events
trackPB_2_new = trackPB_2(1);
for i = 2:length(trackPB_2)
    if meta.data_block1.tvec(trackPB_2(i))- meta.data_block1.tvec(trackPB_2(i-1)) > refractory_period %new event if far away enough in time
        trackPB_2_new(end+1) = trackPB_2(i);
    end
end

trackPB_3 = find(meta.data_block1.channels{2}> threshold); %find possible events
trackPB_3_new = trackPB_3(1);
for i = 2:length(trackPB_3)
    if meta.data_block1.tvec(trackPB_3(i))- meta.data_block1.tvec(trackPB_3(i-1)) > refractory_period %new event if far away enough in time
        trackPB_3_new(end+1) = trackPB_3(i);
    end
end

trackPB_4 = find(meta.data_block1.channels{5}> threshold); %find possible events
trackPB_4_new = trackPB_4(1);
for i = 2:length(trackPB_4)
    if meta.data_block1.tvec(trackPB_4(i))- meta.data_block1.tvec(trackPB_4(i-1)) > refractory_period %new event if far away enough in time
        trackPB_4_new(end+1) = trackPB_4(i);
    end
end

%%
count_PB1 = 2;
trackPB_1_ordered(1) = trackPB_1_new(1);
for i = 2:length(trackPB_1_new)
    if isempty (find(trackPB_2_new > trackPB_1_new(i-1) & trackPB_2_new < trackPB_1_new(i)))
        if isempty (find(trackPB_3_new > trackPB_1_new(i-1) & trackPB_3_new < trackPB_1_new(i)))
            if isempty (find(trackPB_4_new > trackPB_1_new(i-1) & trackPB_4_new < trackPB_1_new(i)))
                continue
            else
                trackPB_1_ordered(count_PB1) = trackPB_1_new(i);
                count_PB1 = count_PB1 + 1;
            end
        else
            trackPB_1_ordered(count_PB1) = trackPB_1_new(i);
            count_PB1 = count_PB1 + 1;
        end
    else
        trackPB_1_ordered(count_PB1) = trackPB_1_new(i);
        count_PB1 = count_PB1 + 1;
    end
end

count_PB2 = 2;
trackPB_2_ordered(1) = trackPB_2_new(1);
for i = 2:length(trackPB_2_new)
    if isempty (find(trackPB_3_new > trackPB_2_new(i-1) & trackPB_3_new < trackPB_2_new(i)))
        if isempty (find(trackPB_4_new > trackPB_2_new(i-1) & trackPB_4_new < trackPB_2_new(i)))
            if isempty (find(trackPB_1_new > trackPB_2_new(i-1) & trackPB_1_new < trackPB_2_new(i)))
                continue
            else
                trackPB_2_ordered(count_PB2) = trackPB_2_new(i);
                count_PB2 = count_PB2 + 1;
            end
        else
            trackPB_2_ordered(count_PB2) = trackPB_2_new(i);
            count_PB2 = count_PB2 + 1;
        end
    else
        trackPB_2_ordered(count_PB2) = trackPB_2_new(i);
        count_PB2 = count_PB2 + 1;
    end
end

count_PB3 = 2;
trackPB_3_ordered(1) = trackPB_3_new(1);
for i = 2:length(trackPB_3_new)
    if isempty (find(trackPB_4_new > trackPB_3_new(i-1) & trackPB_4_new < trackPB_3_new(i)))
        if isempty (find(trackPB_1_new > trackPB_3_new(i-1) & trackPB_1_new < trackPB_3_new(i)))
            if isempty (find(trackPB_2_new > trackPB_3_new(i-1) & trackPB_2_new < trackPB_3_new(i)))
                continue
            else
                trackPB_3_ordered(count_PB3) = trackPB_3_new(i);
                count_PB3 = count_PB3 + 1;
            end
        else
            trackPB_3_ordered(count_PB3) = trackPB_3_new(i);
            count_PB3 = count_PB3 + 1;
        end
    else
        trackPB_3_ordered(count_PB3) = trackPB_3_new(i);
        count_PB3 = count_PB3 + 1;
    end
end

count_PB4 = 2;
trackPB_4_ordered(1) = trackPB_4_new(1);
for i = 2:length(trackPB_4_new)
    if isempty (find(trackPB_1_new > trackPB_4_new(i-1) & trackPB_1_new < trackPB_4_new(i)))
        if isempty (find(trackPB_2_new > trackPB_4_new(i-1) & trackPB_2_new < trackPB_4_new(i)))
            if isempty (find(trackPB_3_new > trackPB_4_new(i-1) & trackPB_3_new < trackPB_4_new(i)))
                continue
            else
                trackPB_4_ordered(count_PB4) = trackPB_4_new(i);
                count_PB4 = count_PB4 + 1;
            end
        else
            trackPB_4_ordered(count_PB4) = trackPB_4_new(i);
            count_PB4 = count_PB4 + 1;
        end
    else
        trackPB_4_ordered(count_PB4) = trackPB_4_new(i);
        count_PB4 = count_PB4 + 1;
    end
end

%%
for ij = 1:4
    switch ij
        case 1
            event_PB1 = find(eventLog_block1.nosepokeID == ij);
        case 2
            event_PB2 = find(eventLog_block1.nosepokeID == ij);
        case 3
            event_PB3 = find(eventLog_block1.nosepokeID == ij);
        case 4
            event_PB4 = find(eventLog_block1.nosepokeID == ij);
    end
end

for kk = 1:length(trackPB_1_ordered)
    if kk > length(event_PB1)
        continue
    else
        offset.PB1(kk) = meta.data_block1.tvec(trackPB_1_ordered(kk)) - eventLog_block1.nosepokeT(event_PB1(kk));
    end
end
for kk = 1:length(trackPB_2_ordered)
    if kk > length(event_PB2)
        continue
    else
        offset.PB2(kk) = meta.data_block1.tvec(trackPB_2_ordered(kk)) - eventLog_block1.nosepokeT(event_PB2(kk));
    end
end
for kk = 1:length(trackPB_3_ordered)
    if kk > length(event_PB3)
        continue
    else
        offset.PB3(kk) = meta.data_block1.tvec(trackPB_3_ordered(kk)) - eventLog_block1.nosepokeT(event_PB3(kk));
    end
end
for kk = 1:length(trackPB_4_ordered)
    if kk > length(event_PB4)
        continue
    else
        offset.PB4(kk) = meta.data_block1.tvec(trackPB_4_ordered(kk)) - eventLog_block1.nosepokeT(event_PB4(kk));
    end
end
for ik = 1:eventLog_block1.trials_initiated
    
    if isempty(find(eventLog_block1.nosepokeT > eventLog_block1.trialT(ik) & eventLog_block1.nosepokeT < eventLog_block1.trialT(ik+1),1,'first'))
        meta.offsetT_block1(ik) = 0;
    else
        offset.temp(ik) = find(eventLog_block1.nosepokeT > eventLog_block1.trialT(ik) & eventLog_block1.nosepokeT < eventLog_block1.trialT(ik+1),1,'first'); %find corresponding photobeam break
        offset.nosepokeID(ik) = eventLog_block1.nosepokeID(offset.temp(ik));
        
        switch offset.nosepokeID(ik)
            case 1
                offset.temp2(ik) = find(event_PB1 == offset.temp(ik),1,'first');
                if offset.temp2(ik) > length(offset.PB1)
                    continue
                else
                    meta.offsetT_block1(ik) = offset.PB1(offset.temp2(ik));
                end
            case 2
                offset.temp2(ik) = find(event_PB2 == offset.temp(ik),1,'first');
                if offset.temp2(ik) > length(offset.PB2)
                    continue
                else
                    meta.offsetT_block1(ik) = offset.PB2(offset.temp2(ik));
                end
            case 3
                offset.temp2(ik) = find(event_PB3 == offset.temp(ik),1,'first');
                if offset.temp2(ik) > length(offset.PB3)
                    continue
                else
                    meta.offsetT_block1(ik) = offset.PB3(offset.temp2(ik));
                end
            case 4
                offset.temp2(ik) = find(event_PB4 == offset.temp(ik),1,'first');
                if offset.temp2(ik) > length(offset.PB4)
                    continue
                else
                    meta.offsetT_block1(ik) = offset.PB4(offset.temp2(ik));
                end
        end
    end
end

block1_startOffset = find(any(meta.offsetT_block1,1),1,'first');
block1_endOffset = find(any(meta.offsetT_block1,1),1,'last');
inc_factor = (meta.offsetT_block1(block1_endOffset) - meta.offsetT_block1(block1_startOffset)) / (block1_endOffset - block1_startOffset);

for jj = 1:length(meta.offsetT_block1)
    if jj == 1
        continue
    end
    if meta.offsetT_block1(jj) == 0
        meta.offsetT_block1(jj) = meta.offsetT_block1(block1_startOffset) + (inc_factor*(jj-1));
    end
end
%% Find time offset between matlab and AMPX  for block 2
threshold = 3000; %threshold for an event
refractory_period = 5; %period of time before passing threshold considered a new event
track2PB_1 = find(meta.data_block2.channels{4}> threshold); %find possible events
track2PB_1_new = track2PB_1(1);
for i = 2:length(track2PB_1)
    if meta.data_block2.tvec(track2PB_1(i))- meta.data_block2.tvec(track2PB_1(i-1)) > refractory_period %new event if far away enough in time
        track2PB_1_new(end+1) = track2PB_1(i);
    end
end

track2PB_2 = find(meta.data_block2.channels{3}> threshold); %find possible events
track2PB_2_new = track2PB_2(1);
for i = 2:length(track2PB_2)
    if meta.data_block2.tvec(track2PB_2(i))- meta.data_block2.tvec(track2PB_2(i-1)) > refractory_period %new event if far away enough in time
        track2PB_2_new(end+1) = track2PB_2(i);
    end
end

track2PB_3 = find(meta.data_block2.channels{2}> threshold); %find possible events
track2PB_3_new = track2PB_3(1);
for i = 2:length(track2PB_3)
    if meta.data_block2.tvec(track2PB_3(i))- meta.data_block2.tvec(track2PB_3(i-1)) > refractory_period %new event if far away enough in time
        track2PB_3_new(end+1) = track2PB_3(i);
    end
end

track2PB_4 = find(meta.data_block2.channels{5}> threshold); %find possible events
track2PB_4_new = track2PB_4(1);
for i = 2:length(track2PB_4)
    if meta.data_block2.tvec(track2PB_4(i))- meta.data_block2.tvec(track2PB_4(i-1)) > refractory_period %new event if far away enough in time
        track2PB_4_new(end+1) = track2PB_4(i);
    end
end

%%
count2_PB1 = 2;
track2PB_1_ordered(1) = track2PB_1_new(1);
for i = 2:length(track2PB_1_new)
    if isempty (find(track2PB_2_new > track2PB_1_new(i-1) & track2PB_2_new < track2PB_1_new(i)))
        if isempty (find(track2PB_3_new > track2PB_1_new(i-1) & track2PB_3_new < track2PB_1_new(i)))
            if isempty (find(track2PB_4_new > track2PB_1_new(i-1) & track2PB_4_new < track2PB_1_new(i)))
                continue
            else
                track2PB_1_ordered(count2_PB1) = track2PB_1_new(i);
                count2_PB1 = count2_PB1 + 1;
            end
        else
            track2PB_1_ordered(count2_PB1) = track2PB_1_new(i);
            count2_PB1 = count2_PB1 + 1;
        end
    else
        track2PB_1_ordered(count2_PB1) = track2PB_1_new(i);
        count2_PB1 = count2_PB1 + 1;
    end
end

count2_PB2 = 2;
track2PB_2_ordered(1) = track2PB_2_new(1);
for i = 2:length(track2PB_2_new)
    if isempty (find(track2PB_3_new > track2PB_2_new(i-1) & track2PB_3_new < track2PB_2_new(i)))
        if isempty (find(track2PB_4_new > track2PB_2_new(i-1) & track2PB_4_new < track2PB_2_new(i)))
            if isempty (find(track2PB_1_new > track2PB_2_new(i-1) & track2PB_1_new < track2PB_2_new(i)))
                continue
            else
                track2PB_2_ordered(count2_PB2) = track2PB_2_new(i);
                count2_PB2 = count2_PB2 + 1;
            end
        else
            track2PB_2_ordered(count2_PB2) = track2PB_2_new(i);
            count2_PB2 = count2_PB2 + 1;
        end
    else
        track2PB_2_ordered(count2_PB2) = track2PB_2_new(i);
        count2_PB2 = count2_PB2 + 1;
    end
end

count2_PB3 = 2;
track2PB_3_ordered(1) = track2PB_3_new(1);
for i = 2:length(track2PB_3_new)
    if isempty (find(track2PB_4_new > track2PB_3_new(i-1) & track2PB_4_new < track2PB_3_new(i)))
        if isempty (find(track2PB_1_new > track2PB_3_new(i-1) & track2PB_1_new < track2PB_3_new(i)))
            if isempty (find(track2PB_2_new > track2PB_3_new(i-1) & track2PB_2_new < track2PB_3_new(i)))
                continue
            else
                track2PB_3_ordered(count2_PB3) = track2PB_3_new(i);
                count2_PB3 = count2_PB3 + 1;
            end
        else
            track2PB_3_ordered(count2_PB3) = track2PB_3_new(i);
            count2_PB3 = count2_PB3 + 1;
        end
    else
        track2PB_3_ordered(count2_PB3) = track2PB_3_new(i);
        count2_PB3 = count2_PB3 + 1;
    end
end

count2_PB4 = 2;
track2PB_4_ordered(1) = track2PB_4_new(1);
for i = 2:length(track2PB_4_new)
    if isempty (find(track2PB_1_new > track2PB_4_new(i-1) & track2PB_1_new < track2PB_4_new(i)))
        if isempty (find(track2PB_2_new > track2PB_4_new(i-1) & track2PB_2_new < track2PB_4_new(i)))
            if isempty (find(track2PB_3_new > track2PB_4_new(i-1) & track2PB_3_new < track2PB_4_new(i)))
                continue
            else
                track2PB_4_ordered(count2_PB4) = track2PB_4_new(i);
                count2_PB4 = count2_PB4 + 1;
            end
        else
            track2PB_4_ordered(count2_PB4) = track2PB_4_new(i);
            count2_PB4 = count2_PB4 + 1;
        end
    else
        track2PB_4_ordered(count2_PB4) = track2PB_4_new(i);
        count2_PB4 = count2_PB4 + 1;
    end
end

%%
for ij = 1:4
    switch ij
        case 1
            event2_PB1 = find(eventLog_block2.nosepokeID == ij);
        case 2
            event2_PB2 = find(eventLog_block2.nosepokeID == ij);
        case 3
            event2_PB3 = find(eventLog_block2.nosepokeID == ij);
        case 4
            event2_PB4 = find(eventLog_block2.nosepokeID == ij);
    end
end

for kk = 1:length(track2PB_1_ordered)
    if kk > length(event2_PB1)
        continue
    else
        offset2.PB1(kk) = meta.data_block2.tvec(track2PB_1_ordered(kk)) - eventLog_block2.nosepokeT(event2_PB1(kk));
    end
end
for kk = 1:length(track2PB_2_ordered)
    if kk > length(event2_PB2)
        continue
    else
        offset2.PB2(kk) = meta.data_block2.tvec(track2PB_2_ordered(kk)) - eventLog_block2.nosepokeT(event2_PB2(kk));
    end
end
for kk = 1:length(track2PB_3_ordered)
    if kk > length(event2_PB3)
        continue
    else
        offset2.PB3(kk) = meta.data_block2.tvec(track2PB_3_ordered(kk)) - eventLog_block2.nosepokeT(event2_PB3(kk));
    end
end
for kk = 1:length(track2PB_4_ordered)
    if kk > length(event2_PB4)
        continue
    else
        offset2.PB4(kk) = meta.data_block2.tvec(track2PB_4_ordered(kk)) - eventLog_block2.nosepokeT(event2_PB4(kk));
    end
end

for ik = 1:eventLog_block2.trials_initiated
    
    if isempty(find(eventLog_block2.nosepokeT > eventLog_block2.trialT(ik) & eventLog_block2.nosepokeT < eventLog_block2.trialT(ik+1),1,'first'))
        meta.offsetT_block2(ik) = 0;
    else
        offset2.temp(ik) = find(eventLog_block2.nosepokeT > eventLog_block2.trialT(ik) & eventLog_block2.nosepokeT < eventLog_block2.trialT(ik+1),1,'first'); %find corresponding photobeam break
        offset2.nosepokeID(ik) = eventLog_block2.nosepokeID(offset2.temp(ik));
        
        switch offset2.nosepokeID(ik)
            case 1
                offset2.temp2(ik) = find(event2_PB1 == offset2.temp(ik),1,'first');
                if offset2.temp2(ik) > length(offset2.PB1)
                    continue
                else
                    meta.offsetT_block2(ik) = offset2.PB1(offset2.temp2(ik));
                end
            case 2
                offset2.temp2(ik) = find(event2_PB2 == offset2.temp(ik),1,'first');
                if offset2.temp2(ik) > length(offset2.PB2)
                    continue
                else
                    meta.offsetT_block2(ik) = offset2.PB2(offset2.temp2(ik));
                end
            case 3
                offset2.temp2(ik) = find(event2_PB3 == offset2.temp(ik),1,'first');
                if offset2.temp2(ik) > length(offset2.PB3)
                    continue
                else
                    meta.offsetT_block2(ik) = offset2.PB3(offset2.temp2(ik));
                end
            case 4
                offset2.temp2(ik) = find(event2_PB4 == offset2.temp(ik),1,'first');
                if offset2.temp2(ik) > length(offset2.PB4)
                    continue
                else
                    meta.offsetT_block2(ik) = offset2.PB4(offset2.temp2(ik));
                end
        end
    end
end

block2_startOffset = find(any(meta.offsetT_block2,1),1,'first');
block2_endOffset = find(any(meta.offsetT_block2,1),1,'last');
inc_factor2 = (meta.offsetT_block2(block2_endOffset) - meta.offsetT_block2(block2_startOffset)) / (block2_endOffset - block2_startOffset);

for jj = 1:length(meta.offsetT_block2)
    if jj == 1
        continue
    end
    if meta.offsetT_block2(jj) == 0
        meta.offsetT_block2(jj) = meta.offsetT_block2(block2_startOffset) + (inc_factor2*(jj-1));
    end
end