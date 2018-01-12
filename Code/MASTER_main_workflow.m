directory = 'E:\Jimmie\Jimmie\Analysis\'; %working directory
destination = 'E:\CueCoding\mat files\'; %where to save .mat files

for ii = 1:4 % work through each rat in the dataset
    switch ii
        case 1
            rat_id = 'R053';
            day = {'-2014-11-04'; '-2014-11-05'; '-2014-11-06'; '-2014-11-07'; '-2014-11-08'; ...
                '-2014-11-10'; '-2014-11-11'; '-2014-11-12'; '-2014-11-13'; '-2014-11-14'; ...
                '-2014-11-15'; '-2014-11-16'};
            block_order_list = [1 1 2 1 2 2 1 2 2 1 1 2];
        case 2
            rat_id = 'R056';
            day = {'-2015-05-16'; '-2015-05-18'; '-2015-05-19'; '-2015-05-21'; '-2015-05-22'; ...
                '-2015-05-23'; '-2015-05-29'; '-2015-05-31'; '-2015-06-01'; '-2015-06-02'; ...
                '-2015-06-03'; '-2015-06-04'; '-2015-06-05'; '-2015-06-08'; '-2015-06-09'};
            block_order_list = [1 1 2 2 2 1 1 2 2 1 2 1 2 1 2];
        case 3
            rat_id = 'R057';
            day = {'-2015-02-14'; '-2015-02-15'; '-2015-02-16'; '-2015-02-17'; '-2015-02-18'; ...
                '-2015-02-21'; '-2015-02-22'; '-2015-02-23'; '-2015-02-24'; '-2015-02-25'; ...
                '-2015-02-26'; '-2015-02-27'};
            block_order_list = [2 1 2 1 2 1 2 1 2 1 2 1];
        case 4
            rat_id = 'R060';
            day = {'-2014-12-11'; '-2014-12-12'; '-2014-12-13'; '-2014-12-14'; '-2014-12-15'; ...
                '-2014-12-17'; '-2014-12-19'; '-2014-12-20'; '-2014-12-21'; '-2014-12-23'; '-2014-12-24'; ...
                '-2014-12-26'; '-2014-12-27'; '-2014-12-28'; '-2014-12-29'; '-2014-12-30'; '-2014-12-31'; ...
                '-2015-01-01'; '-2015-01-02'; '-2015-01-03'; '-2015-01-04'};
            block_order_list = [1 2 1 2 1 2 1 2 1 1 2 1 1 2 1 1 2 1 2 1 2];
    end
    for kj = 1:length(day)
        
        %% input session info
        sesh.session_id = cat(2,rat_id,day{kj}); % current day for analysis
        
        sesh.block_order = block_order_list(kj); % 1 if light block came first, 2 if sound block came first
        disp(cat(2,'loading session ',sesh.session_id));
        switch ii
            case 1
                meta = metaFnc_R053(sesh.session_id,sesh.block_order); % for R053
            case 2
                meta = metaFnc_new_RR2(sesh.session_id,sesh.block_order); % for R056
            otherwise
                meta = metaFnc_old_RR2(sesh.session_id,sesh.block_order); % for R057 and R060
        end
        
        metadata.TrialInfo{1} = meta.TrialInfo_block1;
        metadata.TrialInfo{2} = meta.TrialInfo_block2;
        metadata.TrialInfo{1}.offsetT = meta.offsetT_block1;
        metadata.TrialInfo{2}.offsetT = meta.offsetT_block2;
        
        for kk = 1:16
            %% input tt info
            sesh.tt_number = kk; % which tt are you currently looking at
            disp(cat(2,'tetrode ',num2str(kk)));
            cd(cat(2,directory,sesh.session_id(1:4),'\',sesh.session_id,'\ntt\TT',num2str(sesh.tt_number)));
            if isempty(dir('*.t')) == 1
                continue
            end
            
            for jj = 1:length(dir('*.t'))
                t_files = dir('*.t');
                sesh.cell_number = jj; % which cell number on current tt
                sesh.tt_fname = t_files(jj).name;
                
                %% load spikes
                disp('loading spikes')
                spk_t = LoadSpikes({sesh.tt_fname}); % dbl check 32 v 64 bit
                spk_t = Data(spk_t{1});
                
                %% get event times
                dataPoint = genDP(meta);
                metadata.dataPoint = dataPoint;
                
                %% peri-event histograms
                sesh.PETH.Trial = 1; %generate PETH separated by cue conditions at trial start
                sesh.PETH.Arm = 1; %generate PETH separated by arm location
                sesh.PETH.Approach = 1; %generate PETH separated by approach and skip for unrewarded trials
                sesh.PETH.Trial_np = 1; %generate PETH separated by cue conditions at trial start for when the rat approached
                
                PETH = genPETH(sesh,meta,spk_t,dataPoint);
                
                %% rasterplots
                RAST = genRAST(meta,spk_t,dataPoint);
                
                %% extract event related firing rates
                FRATE = genFRATE(meta,spk_t,dataPoint);

                %% bootstrap
                cfg_in.method = 2;
                SHUFF.Trial = genSHUFF(cfg_in,spk_t,FRATE.Interval.Trial);
                
                %% diff
                trial_int_length = (FRATE.Interval.Trial.tend - FRATE.Interval.Trial.tstart)/2;
                trial_event_times = dataPoint.Trials/1000;
                
                cfg_in = [];
                DIFF.Trial = genDIFF(cfg_in,SHUFF.Trial,trial_event_times,trial_int_length);
                
                %% rank
                RANK.Trial = genRANK(DIFF.Trial,FRATE.Task.Trial_firing_rate,FRATE.Task.Trial_B4_firing_rate);
                
                %% Wilcoxon signed rank and Mann-whitney u-tests
                TESTS = genTESTS(FRATE);
                
                %% save variables in .mat file
                save(cat(2,destination,sesh.session_id,'-TT',num2str(sesh.tt_number),'-cell',num2str(sesh.cell_number)),'metadata','FRATE','TESTS','PETH','RAST','spk_t','dataPoint','sesh','SHUFF','DIFF','RANK');
                clearvars -except ii kj kk jj sesh rat_id days block_order_list directory destination meta
            end
        end
    end
end