%% R060 light

data_days = 70; % how many days of data
cum_rolling_avg_light = linspace(0,0,200);
first_vs_second_half = 1;

for ij = 1:data_days
    clear behav_analysis
    switch ij
        case 1
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-16\R060-2014-07-16_light_behaviour.mat')
case 2
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-19\R060-2014-07-19_light_behaviour.mat')
case 3
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-21\R060-2014-07-21_light_behaviour.mat')
case 4
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-22\R060-2014-07-22_light_behaviour.mat')
case 5
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-23\R060-2014-07-23_light_behaviour.mat')
case 6
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-24\R060-2014-07-24_light_behaviour.mat')
case 7
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-25\R060-2014-07-25_light_behaviour.mat')
case 8
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-26\R060-2014-07-26_light_behaviour.mat')
case 9
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-27\R060-2014-07-27_light_behaviour.mat')
case 10
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-30\R060-2014-07-30_light_behaviour.mat')
case 11
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-31\R060-2014-07-31_light_behaviour.mat')
case 12
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-01\R060-2014-08-01_light_behaviour.mat')
case 13
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-02\R060-2014-08-02_light_behaviour.mat')
case 14
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-03\R060-2014-08-03_light_behaviour.mat')
case 15
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-04\R060-2014-08-04_light_behaviour.mat')
case 16
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-05\R060-2014-08-05_light_behaviour.mat')
case 17
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-06\R060-2014-08-06_light_behaviour.mat')
        case 18
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-08\R060-2014-08-08_light_behaviour.mat')
        case 19
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-09\R060-2014-08-09_light_behaviour.mat')
        case 20
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-11\R060-2014-08-11_light_behaviour.mat')
        case 21
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-12\R060-2014-08-12_light_behaviour.mat')
%         case 22 no data
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-13\R060-2014-08-13_light_behaviour.mat')
        case 22
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-14\R060-2014-08-14_light_behaviour.mat')
%         case 24 no data
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-15\R060-2014-08-15_light_behaviour.mat')
        case 23
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-16\R060-2014-08-16_light_behaviour.mat')
        case 24
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-18\R060-2014-08-18_light_behaviour.mat')
%         case 26 - data didn't save (no data)
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-19\R060-2014-08-19_light_behaviour.mat')
        case 25
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-20\R060-2014-08-20_light_behaviour.mat')
        case 26
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-21\R060-2014-08-21_light_behaviour.mat')
        case 27
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-22\R060-2014-08-22_light_behaviour.mat')
        case 28
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-23\R060-2014-08-23_light_behaviour.mat')
        case 29
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-24\R060-2014-08-24_light_behaviour.mat')
        case 30
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-25\R060-2014-08-25_light_behaviour.mat')
        case 31
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-26\R060-2014-08-26_light_behaviour.mat')
        case 32
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-27\R060-2014-08-27_light_behaviour.mat')
        case 33
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-28\R060-2014-08-28_light_behaviour.mat')
        case 34
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-29\R060-2014-08-29_light_behaviour.mat')
        case 35
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-30\R060-2014-08-30_light_behaviour.mat')
        case 36
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-31\R060-2014-08-31_light_behaviour.mat')
        case 37
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-01\R060-2014-09-01_light_behaviour.mat')
        case 38
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-02\R060-2014-09-02_light_behaviour.mat')
        case 39
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-03\R060-2014-09-03_light_behaviour.mat')
        case 40
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-05\R060-2014-09-05_light_behaviour.mat')
            case 41
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-06\R060-2014-09-06_light_behaviour.mat')
case 42 % recording
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-11\R060-2014-12-11-light_behaviour.mat')
        case 43
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-12\R060-2014-12-12-light_behaviour.mat')
        case 44
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-13\R060-2014-12-13-light_behaviour.mat')
        case 45
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-14\R060-2014-12-14-light_behaviour.mat')
        case 46
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-15\R060-2014-12-15-light_behaviour.mat')
        case 47
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-17\R060-2014-12-17-light_behaviour.mat')
        case 48
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-19\R060-2014-12-19-light_behaviour.mat')
        case 49
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-20\R060-2014-12-20-light_behaviour.mat')
        case 50
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-21\R060-2014-12-21-light_behaviour.mat')
        case 51
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-22\R060-2014-12-22-light_behaviour.mat')
        case 52
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-23\R060-2014-12-23-light_behaviour.mat')
        case 53
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-24\R060-2014-12-24-light_behaviour.mat')
        case 54
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-26\R060-2014-12-26-light_behaviour.mat')
        case 55
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-27\R060-2014-12-27-light_behaviour.mat')
        case 56
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-28\R060-2014-12-28-light_behaviour.mat')
        case 57
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-29\R060-2014-12-29-light_behaviour.mat')
        case 58
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-30\R060-2014-12-30-light_behaviour.mat')
        case 59
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-31\R060-2014-12-31-light_behaviour.mat')
        case 60
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-01\R060-2015-01-01-light_behaviour.mat')
        case 61
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-02\R060-2015-01-02-light_behaviour.mat')
        case 62
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-03\R060-2015-01-03-light_behaviour.mat')
        case 63
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-04\R060-2015-01-04-light_behaviour.mat')
        case 64
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-05\R060-2015-01-05-light_behaviour.mat')
        case 65
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-06\R060-2015-01-06-light_behaviour.mat')
        case 66
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-07\R060-2015-01-07-light_behaviour.mat')
        case 67
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-08\R060-2015-01-08-light_behaviour.mat')
        case 68
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-09\R060-2015-01-09-light_behaviour.mat')
        case 69
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-10\R060-2015-01-10-light_behaviour.mat')
        case 70
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-11\R060-2015-01-11-light_behaviour.mat')
%         case 71 % both block introduced
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-12\R060-2015-01-12-light_behaviour.mat')
%         case 72
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-13\R060-2015-01-13-light_behaviour.mat')
%         case 73
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-14\R060-2015-01-14-light_behaviour.mat')
%         case 74
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-15\R060-2015-01-15-light_behaviour.mat')
%         case 75
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-16\R060-2015-01-16-light_behaviour.mat')
%         case 76
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-17\R060-2015-01-17-light_behaviour.mat')
%         case 77
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-18\R060-2015-01-18-light_behaviour.mat')
%         case 78
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-19\R060-2015-01-19-light_behaviour.mat')
%         case 79
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-20\R060-2015-01-20-light_behaviour.mat')
%         case 80
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-21\R060-2015-01-21-light_behaviour.mat')
%         case 81
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-22\R060-2015-01-22-light_behaviour.mat')

    end
    
    behav_analysis.total_approached = 0;
    behav_analysis.total_skipped = 0;
    behav_analysis.trials = eventLog.trials_initiated;
    nTrials(ij) = eventLog.trials_initiated;
    behav_analysis.light_1 = [];
    behav_analysis.light_2 = [];
    ind = find(isnan(eventLog.trialT));
    eventLog.trialT(ind) = [];
    
    for i = 1:length(eventLog.trialT) %length of trials
        temp = find (eventLog.pb_breaksT == eventLog.trialT(i),1,'first'); %find the corresponding time in eventLog.pb_breaksT of the trial initiation
        if eventLog.pb_breaksID(temp + 1) == eventLog.pb_breaksID(temp) + 4 %see if the next pb break is the adjacent reward pb
            temp2 = find (eventLog.ALLnosepokeT == eventLog.pb_breaksT(temp +1),1,'first'); %find the location of this nosepoke in the eventLog.ALLnosepokeT variable
            for jj = temp2:length(eventLog.unnosepokeT)
                if eventLog.ALLnosepokeID(temp2) ~= eventLog.unnosepokeID(jj) %find the last unnosepoke from current reward pb
                    temp3 = eventLog.unnosepokeT(jj - 1) - eventLog.ALLnosepokeT(temp2); %find time difference from first nosepoke to last unnosepoke
                    break %your work here is done, move on
                end
            end
            if i > 100 % old part of code that manages variable that separates the task into two halves
                first_vs_second_half = 1;
            else
                first_vs_second_half = 0;
            end
            if temp3 > 1
                behav_analysis.approached(i) = 1;
                behav_analysis.total_approached = behav_analysis.total_approached + 1;
                behav_analysis.approached_ID(behav_analysis.total_approached) = eventLog.pb_breaksID(temp);
                behav_analysis.light_approached(behav_analysis.total_approached) = eventLog.light_ID(i);
               
                switch eventLog.light_ID(i)
                    case 1
                        behav_analysis.light_1(length(behav_analysis.light_1) + 1) = 1;
                        behav_analysis.light_1_half(length(behav_analysis.light_1)) = first_vs_second_half;
                    case 2
                        behav_analysis.light_2(length(behav_analysis.light_2) + 1) = 1;
                        behav_analysis.light_2_half(length(behav_analysis.light_2)) = first_vs_second_half;
                end             
                
            else
                behav_analysis.approached(i) = 0;
                behav_analysis.total_skipped = behav_analysis.total_skipped + 1;
                behav_analysis.skipped_ID(behav_analysis.total_skipped) = eventLog.pb_breaksID(temp);
                behav_analysis.light_skipped(behav_analysis.total_skipped) = eventLog.light_ID(i);
                switch eventLog.light_ID(i)
                    case 1
                        behav_analysis.light_1(length(behav_analysis.light_1) + 1) = 0;
                        behav_analysis.light_1_half(length(behav_analysis.light_1)) = first_vs_second_half;
                    case 2
                        behav_analysis.light_2(length(behav_analysis.light_2) + 1) = 0;
                        behav_analysis.light_2_half(length(behav_analysis.light_2)) = first_vs_second_half;
                end
                
            end
        else
            behav_analysis.approached(i) = 0;
            behav_analysis.total_skipped = behav_analysis.total_skipped + 1;
            behav_analysis.skipped_ID(behav_analysis.total_skipped) = eventLog.pb_breaksID(temp);
            behav_analysis.light_skipped(behav_analysis.total_skipped) = eventLog.light_ID(i);
            switch eventLog.light_ID(i)
                case 1
                    behav_analysis.light_1(length(behav_analysis.light_1) + 1) = 0;
                    behav_analysis.light_1_half(length(behav_analysis.light_1)) = first_vs_second_half;
                case 2
                    behav_analysis.light_2(length(behav_analysis.light_2) + 1) = 0;
                    behav_analysis.light_2_half(length(behav_analysis.light_2)) = first_vs_second_half;
            end

        end
        
        behav_analysis.count_approached(i) = behav_analysis.total_approached;
        behav_analysis.count_skipped(i) = behav_analysis.total_skipped;
        behav_analysis.percent_skipped = (behav_analysis.total_skipped ./ i) .* 100;
        behav_analysis.rewarded(i) = eventLog.rewarded(i);
        
    end
    

    for j = 1:2
        temp = find(behav_analysis.light_approached == j);
        behav_analysis.total_light_approached(j) = length(temp);
        if behav_analysis.total_skipped == 0;
        else
        temp = find(behav_analysis.light_skipped == j);
        behav_analysis.total_light_skipped(j) = length(temp);
        end

    end
    
    light_app_1(ij) = behav_analysis.total_light_approached(1);
    light_app_2(ij) = behav_analysis.total_light_approached(2);
    if behav_analysis.total_skipped == 0;
        light_skip_1(ij) = 0;
        light_skip_2(ij) = 0;
    else
    light_skip_1(ij) = behav_analysis.total_light_skipped(1);
    light_skip_2(ij) = behav_analysis.total_light_skipped(2);
    end
    
    behav_analysis.rolling_avg = linspace(0,0,200);
    
    for k = 1:length(eventLog.trialT)
        if behav_analysis.rewarded(k) == 2
            behav_analysis.rewarded(k) = 1;
        end
        if behav_analysis.approached(k) == behav_analysis.rewarded(k)
            behav_analysis.correct(k) = 1;
        else
            behav_analysis.correct(k) = 0;
        end
        if k > 4
            behav_analysis.rolling_avg(k-4) = mean(behav_analysis.correct(k-4:k));
        end
    end
        
    cum_rolling_avg_light = cum_rolling_avg_light + behav_analysis.rolling_avg;
    prop_correct(ij) = mean(behav_analysis.correct);
end

prop_light_app_1 = light_app_1 ./ (light_app_1 + light_skip_1);
prop_light_app_2 = light_app_2 ./ (light_app_2 + light_skip_2);
prop_light_skip_1 = light_skip_1 ./ (light_app_1 + light_skip_1);
prop_light_skip_2 = light_skip_2 ./ (light_app_2 + light_skip_2);
rolling_avg_light = cum_rolling_avg_light / ij;

% R053
% chi square test
for i_day = 1:length(light_app_1)
    total_nosepokes(i_day) = light_app_1(i_day) + light_app_2(i_day);
    total_light_1_trials(i_day) = light_app_1(i_day) + light_skip_1(i_day);
    total_light_2_trials(i_day) = light_app_2(i_day) + light_skip_2(i_day);
    expected_pokes_light_1(i_day) = (total_light_1_trials(i_day) * total_nosepokes(i_day)) / nTrials(i_day);
    expected_pokes_light_2(i_day) = (total_light_2_trials(i_day) * total_nosepokes(i_day)) / nTrials(i_day);
    chi_light_1(i_day) = (light_app_1(i_day) - expected_pokes_light_1(i_day))^2 / expected_pokes_light_1(i_day);
    chi_light_2(i_day) = (light_app_2(i_day) - expected_pokes_light_2(i_day))^2 / expected_pokes_light_2(i_day);
    p_light_1(i_day) = 1 - chi2cdf(chi_light_1(i_day),1);
    p_light_2(i_day) = 1 - chi2cdf(chi_light_2(i_day),1);
    
end

%% R060 sound
% clear
data_days = 70; % how many days of data
cum_rolling_avg_sound = linspace(0,0,200);
first_vs_second_half = 1;

for ij = 1:data_days
    clear behav_analysis
    switch ij
        case 1
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-16\R060-2014-07-16_sound_behaviour.mat')
case 2
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-19\R060-2014-07-19_sound_behaviour.mat')
case 3
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-21\R060-2014-07-21_sound_behaviour.mat')
case 4
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-22\R060-2014-07-22_sound_behaviour.mat')
case 5
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-23\R060-2014-07-23_sound_behaviour.mat')
case 6
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-24\R060-2014-07-24_sound_behaviour.mat')
case 7
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-25\R060-2014-07-25_sound_behaviour.mat')
case 8
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-26\R060-2014-07-26_sound_behaviour.mat')
case 9
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-27\R060-2014-07-27_sound_behaviour.mat')
case 10
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-30\R060-2014-07-30_sound_behaviour.mat')
case 11
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-07-31\R060-2014-07-31_sound_behaviour.mat')
% case 12
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-01\R060-2014-08-01_sound_behaviour.mat')
case 12
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-02\R060-2014-08-02_sound_behaviour.mat')
case 13
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-03\R060-2014-08-03_sound_behaviour.mat')
case 14
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-04\R060-2014-08-04_sound_behaviour.mat')
case 15
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-05\R060-2014-08-05_sound_behaviour.mat')
case 16
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-06\R060-2014-08-06_sound_behaviour.mat')
        case 17
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-08\R060-2014-08-08_sound_behaviour.mat')
%         case 18 - data didn't save
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-09\R060-2014-08-09_sound_behaviour.mat')
        case 18
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-11\R060-2014-08-11_sound_behaviour.mat')
        case 19
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-12\R060-2014-08-12_sound_behaviour.mat')
        case 20
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-13\R060-2014-08-13_sound_behaviour.mat')
        case 21
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-14\R060-2014-08-14_sound_behaviour.mat')
        case 22
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-15\R060-2014-08-15_sound_behaviour.mat')
        case 23
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-16\R060-2014-08-16_sound_behaviour.mat')
       case 24
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-18\R060-2014-08-18_sound_behaviour.mat')
        case 25
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-19\R060-2014-08-19_sound_behaviour.mat')
        case 26
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-20\R060-2014-08-20_sound_behaviour.mat')
        case 27
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-21\R060-2014-08-21_sound_behaviour.mat')
        case 28
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-22\R060-2014-08-22_sound_behaviour.mat')
        case 29
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-23\R060-2014-08-23_sound_behaviour.mat')
        case 30
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-24\R060-2014-08-24_sound_behaviour.mat')
        case 31
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-25\R060-2014-08-25_sound_behaviour.mat')
        case 32
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-26\R060-2014-08-26_sound_behaviour.mat')
        case 33
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-27\R060-2014-08-27_sound_behaviour.mat')
        case 34
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-28\R060-2014-08-28_sound_behaviour.mat')
        case 35
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-29\R060-2014-08-29_sound_behaviour.mat')
        case 36
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-30\R060-2014-08-30_sound_behaviour.mat')
        case 37
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-08-31\R060-2014-08-31_sound_behaviour.mat')
        case 38
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-01\R060-2014-09-01_sound_behaviour.mat')
        case 39
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-02\R060-2014-09-02_sound_behaviour.mat')
        case 40
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-03\R060-2014-09-03_sound_behaviour.mat')
        case 41
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-05\R060-2014-09-05_sound_behaviour.mat')
            case 42
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-09-06\R060-2014-09-06_sound_behaviour.mat')
 case 43 % recording
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-11\R060-2014-12-11-sound_behaviour.mat')
        case 44
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-12\R060-2014-12-12-sound_behaviour.mat')
        case 45
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-13\R060-2014-12-13-sound_behaviour.mat')
        case 46
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-14\R060-2014-12-14-sound_behaviour.mat')
        case 47
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-15\R060-2014-12-15-sound_behaviour.mat')
        case 48
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-17\R060-2014-12-17-sound_behaviour.mat')
        case 49
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-19\R060-2014-12-19-sound_behaviour.mat')
        case 50
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-20\R060-2014-12-20-sound_behaviour.mat')
        case 51
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-21\R060-2014-12-21-sound_behaviour.mat')
%         case 10
%              load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-22\R060-2014-12-22-sound_behaviour.mat')
        case 52
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-23\R060-2014-12-23-sound_behaviour.mat')
        case 53
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-24\R060-2014-12-24-sound_behaviour.mat')
        case 54
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-26\R060-2014-12-26-sound_behaviour.mat')
        case 55
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-27\R060-2014-12-27-sound_behaviour.mat')
        case 56
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-28\R060-2014-12-28-sound_behaviour.mat')
        case 57
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-29\R060-2014-12-29-sound_behaviour.mat')
        case 58
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-30\R060-2014-12-30-sound_behaviour.mat')
        case 59
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2014-12-31\R060-2014-12-31-sound_behaviour.mat')
        case 60
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-01\R060-2015-01-01-sound_behaviour.mat')
        case 61
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-02\R060-2015-01-02-sound_behaviour.mat')
        case 62
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-03\R060-2015-01-03-sound_behaviour.mat')
        case 63
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-04\R060-2015-01-04-sound_behaviour.mat')
        case 64
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-05\R060-2015-01-05-sound_behaviour.mat')
        case 65
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-06\R060-2015-01-06-sound_behaviour.mat')
        case 66
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-07\R060-2015-01-07-sound_behaviour.mat')
        case 67
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-08\R060-2015-01-08-sound_behaviour.mat')
        case 68
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-09\R060-2015-01-09-sound_behaviour.mat')
        case 69
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-10\R060-2015-01-10-sound_behaviour.mat')
        case 70
            load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-11\R060-2015-01-11-sound_behaviour.mat')
%         case 71 % both block introduced
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-12\R060-2015-01-12-sound_behaviour.mat')
%         case 72
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-13\R060-2015-01-13-sound_behaviour.mat')
%         case 73
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-14\R060-2015-01-14-sound_behaviour.mat')
%         case 74
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-15\R060-2015-01-15-sound_behaviour.mat')
%         case 75
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-16\R060-2015-01-16-sound_behaviour.mat')
%         case 76
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-17\R060-2015-01-17-sound_behaviour.mat')
%         case 77
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-18\R060-2015-01-18-sound_behaviour.mat')
%         case 78
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-19\R060-2015-01-19-sound_behaviour.mat')
%         case 79
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-20\R060-2015-01-20-sound_behaviour.mat')
%         case 80
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-21\R060-2015-01-21-sound_behaviour.mat')
%         case 81
%             load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-01-22\R060-2015-01-22-sound_behaviour.mat')
%      
    end
    
    behav_analysis.total_approached = 0;
    behav_analysis.total_skipped = 0;
    behav_analysis.trials = eventLog.trials_initiated;
    nTrials(ij) = eventLog.trials_initiated;
    behav_analysis.sound_1 = [];
    behav_analysis.sound_2 = [];
    ind = find(isnan(eventLog.trialT));
    eventLog.trialT(ind) = [];
    
    for i = 1:length(eventLog.trialT)
        temp = find (eventLog.pb_breaksT == eventLog.trialT(i),1,'first');
        if eventLog.pb_breaksID(temp + 1) == eventLog.pb_breaksID(temp) + 4
            temp2 = find (eventLog.ALLnosepokeT == eventLog.pb_breaksT(temp +1),1,'first');
            for jj = temp2:length(eventLog.unnosepokeT)
                if eventLog.ALLnosepokeID(temp2) ~= eventLog.unnosepokeID(jj)
                    temp3 = eventLog.unnosepokeT(jj - 1) - eventLog.ALLnosepokeT(temp2);
                    break
                end
            end
            if i > 100
                first_vs_second_half = 1;
            else
                first_vs_second_half = 0;
            end
            if temp3 > 1
                behav_analysis.approached(i) = 1;
                behav_analysis.total_approached = behav_analysis.total_approached + 1;
                behav_analysis.approached_ID(behav_analysis.total_approached) = eventLog.pb_breaksID(temp);
                behav_analysis.sound_approached(behav_analysis.total_approached) = eventLog.sound_ID(i);
                switch eventLog.sound_ID(i)
                    case 1
                        behav_analysis.sound_1(length(behav_analysis.sound_1) + 1) = 1;
                        behav_analysis.sound_1_half(length(behav_analysis.sound_1)) = first_vs_second_half;
                    case 2
                        behav_analysis.sound_2(length(behav_analysis.sound_2) + 1) = 1;
                        behav_analysis.sound_2_half(length(behav_analysis.sound_2)) = first_vs_second_half;
                end
                
            else
                behav_analysis.approached(i) = 0;
                behav_analysis.total_skipped = behav_analysis.total_skipped + 1;
                behav_analysis.skipped_ID(behav_analysis.total_skipped) = eventLog.pb_breaksID(temp);
                behav_analysis.sound_skipped(behav_analysis.total_skipped) = eventLog.sound_ID(i);
                switch eventLog.sound_ID(i)
                    case 1
                        behav_analysis.sound_1(length(behav_analysis.sound_1) + 1) = 0;
                        behav_analysis.sound_1_half(length(behav_analysis.sound_1)) = first_vs_second_half;
                    case 2
                        behav_analysis.sound_2(length(behav_analysis.sound_2) + 1) = 0;
                        behav_analysis.sound_2_half(length(behav_analysis.sound_2)) = first_vs_second_half;
                end
            end
        else
            behav_analysis.approached(i) = 0;
            behav_analysis.total_skipped = behav_analysis.total_skipped + 1;
            behav_analysis.skipped_ID(behav_analysis.total_skipped) = eventLog.pb_breaksID(temp);
            behav_analysis.sound_skipped(behav_analysis.total_skipped) = eventLog.sound_ID(i);
            switch eventLog.sound_ID(i)
                case 1
                    behav_analysis.sound_1(length(behav_analysis.sound_1) + 1) = 0;
                    behav_analysis.sound_1_half(length(behav_analysis.sound_1)) = first_vs_second_half;
                case 2
                    behav_analysis.sound_2(length(behav_analysis.sound_2) + 1) = 0;
                    behav_analysis.sound_2_half(length(behav_analysis.sound_2)) = first_vs_second_half;
            end
        end
        
        behav_analysis.count_approached(i) = behav_analysis.total_approached;
        behav_analysis.count_skipped(i) = behav_analysis.total_skipped;
        behav_analysis.percent_skipped = (behav_analysis.total_skipped ./ i) .* 100;
        behav_analysis.rewarded(i) = eventLog.rewarded(i);
    end
    
    for j = 1:2
        temp = find(behav_analysis.sound_approached == j);
        behav_analysis.total_sound_approached(j) = length(temp);
        if behav_analysis.total_skipped == 0;
        else
        temp = find(behav_analysis.sound_skipped == j);
        behav_analysis.total_sound_skipped(j) = length(temp);
        end
    end
    
    sound_app_1(ij) = behav_analysis.total_sound_approached(1);
    sound_app_2(ij) = behav_analysis.total_sound_approached(2);
        if behav_analysis.total_skipped == 0;
        sound_skip_1(ij) = 0;
        sound_skip_2(ij) = 0;
    else
    sound_skip_1(ij) = behav_analysis.total_sound_skipped(1);
    sound_skip_2(ij) = behav_analysis.total_sound_skipped(2);
        end
    behav_analysis.rolling_avg = linspace(0,0,200);
    
    for k = 1:length(eventLog.trialT)
        if behav_analysis.rewarded(k) == 2
            behav_analysis.rewarded(k) = 1;
        end
        if behav_analysis.approached(k) == behav_analysis.rewarded(k)
            behav_analysis.correct(k) = 1;
        else
            behav_analysis.correct(k) = 0;
        end
        if k > 4
            behav_analysis.rolling_avg(k-4) = mean(behav_analysis.correct(k-4:k));
        end
    end
        
    cum_rolling_avg_sound = cum_rolling_avg_sound + behav_analysis.rolling_avg;
     prop_correct(ij) = mean(behav_analysis.correct);
end

prop_sound_app_1 = sound_app_1 ./ (sound_app_1 + sound_skip_1);
prop_sound_app_2 = sound_app_2 ./ (sound_app_2 + sound_skip_2);
prop_sound_skip_1 = sound_skip_1 ./ (sound_app_1 + sound_skip_1);
prop_sound_skip_2 = sound_skip_2 ./ (sound_app_2 + sound_skip_2);
rolling_avg_sound = cum_rolling_avg_sound / ij;

% chi square test
for i_day = 1:length(sound_app_1)
    total_nosepokes(i_day) = sound_app_1(i_day) + sound_app_2(i_day);
    total_sound_1_trials(i_day) = sound_app_1(i_day) + sound_skip_1(i_day);
    total_sound_2_trials(i_day) = sound_app_2(i_day) + sound_skip_2(i_day);
    expected_pokes_sound_1(i_day) = (total_sound_1_trials(i_day) * total_nosepokes(i_day)) / nTrials(i_day);
    expected_pokes_sound_2(i_day) = (total_sound_2_trials(i_day) * total_nosepokes(i_day)) / nTrials(i_day);
    chi_sound_1(i_day) = (sound_app_1(i_day) - expected_pokes_sound_1(i_day))^2 / expected_pokes_sound_1(i_day);
    chi_sound_2(i_day) = (sound_app_2(i_day) - expected_pokes_sound_2(i_day))^2 / expected_pokes_sound_2(i_day);
    p_sound_1(i_day) = 1 - chi2cdf(chi_sound_1(i_day),1);
    p_sound_2(i_day) = 1 - chi2cdf(chi_sound_2(i_day),1);
end

%%
% figure
% plot(prop_light_app_1,'color','r'); title('Light'); ylabel('Proportion approached');
% ylim([0 1.1]); hold on; plot(prop_light_app_2,'color','g'); 
% plot(42,0:.05:1,'.','color','black'); xlim([0 70]); xlabel('Session number');
% box off;
% 
% light_sig_first = find(p_light_1 < .05, 1,'first')
% light_sig_first_last = find(p_light_1(light_sig_first + 1:end) > .05, 1,'first') + light_sig_first
% light_sig_next = find(p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_first_last
% light_sig_next_last = find(p_light_1(light_sig_next + 1:end) > .05, 1,'first') + light_sig_next
% light_sig_final = find(p_light_1(light_sig_first_last + 1:end) < .05, 1,'first') + light_sig_next_last
% light_sig_final_last = length(p_light_1)
% plot([light_sig_first light_sig_first_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% plot([light_sig_next light_sig_next_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% plot([light_sig_final light_sig_final_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% 
% figure
% plot(prop_sound_app_1,'color','c'); title('Sound'); ylim([0 1.1]); %xlabel('Session number'); 
% hold on; plot(prop_sound_app_2,'color','b'); plot(42,0:.05:1,'.','color','black'); 
% xlim([0 70]); xlabel('Session number'); ylabel('Proportion approached'); 
% box off;
% 
% sound_sig_first = find(p_sound_1 < .05, 1,'first')
% sound_sig_first_last = find(p_sound_1(sound_sig_first + 1:end) > .05, 1,'first') + sound_sig_first
% sound_sig_next = find(p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_first_last
% sound_sig_next_last = find(p_sound_1(sound_sig_next + 1:end) > .05, 1,'first') + sound_sig_next
% sound_sig_final = find(p_sound_1(sound_sig_first_last + 1:end) < .05, 1,'first') + sound_sig_next_last
% sound_sig_final_last = length(p_sound_1)
% plot([sound_sig_first sound_sig_first_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% plot([sound_sig_next sound_sig_next_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% plot([sound_sig_final sound_sig_final_last], [1.025 1.025], '-k', 'LineWidth',2,'color','r')
% 
% % figure(11)
% % subplot(2,1,1); plot(p_light_1); title('Light rewarded'); ylabel('p value'); 
% % hold on; plot(1:1:length(light_app_1),.05,'.','color',[1 0 0]); hold on; plot(42,0:.05:1,'.','color','black'); 
% % 
% % subplot(2,1,2); plot(p_sound_2); title('Tone rewarded'); ylabel('p value'); xlabel('session');
% % hold on; plot(1:1:length(sound_app_2),.05,'.','color',[1 0 0]);  hold on; plot(42,0:.05:1,'.','color','black'); 
