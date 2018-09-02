function PROCESS = genPROCESS(directory,destination,PETH_generation)
% function PROCESS = genPROCESS(directory,destination,PETH_generation)
%
%
% INPUTS:
%
% OUTPUTS:
%
f = filesep

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
            day = {'-2014-12-12'; '-2014-12-13'; '-2014-12-14'; '-2014-12-15'; ...
                '-2014-12-17'; '-2014-12-20'; '-2014-12-23'; '-2014-12-24'; ...
                '-2014-12-26'; '-2014-12-27'; '-2014-12-28'; '-2014-12-29'; '-2014-12-30'; '-2014-12-31'; ...
                '-2015-01-01'; '-2015-01-02'; '-2015-01-03'; '-2015-01-04'};
            block_order_list = [1 2 1 2 1 2 1 2 1 1 2 1 1 2 1 1 2 1 2 1 2];
    end
    for kj = 1:length(day)
        
        %% input session info
        sesh.session_id = cat(2,rat_id,day{kj}); % current day for analysis
        
        sesh.block_order = block_order_list(kj); % 1 if light block came first, 2 if sound block came first
        disp(cat(2,'loading session ',sesh.session_id));
        cd(cat(2,directory,sesh.session_id(1:4),f,sesh.session_id))
        load(strcat(sesh.session_id,'_metadata.mat'))
        
        %% find .t files
        t_files = dir('*.t');
        for kk = 1:16
            if kk > 9
                str_length = 21;
            else
                str_length = 20;
            end
            %% input tt info
            sesh.tt_number = kk; % which tt are you currently looking at
            disp(cat(2,'tetrode ',num2str(kk)));
            if isempty(dir(strcat(sesh.session_id,'-TT',num2str(kk),'*'))) == 1
                continue
            end
            
            cell_count = 1;
            
            for jj = 1:length(dir('*.t'))
                if strcmp(t_files(jj).name(1:str_length),strcat(sesh.session_id,'-TT',num2str(kk),'_')) == 1
                    sesh.cell_number = cell_count; % which cell number on current tt
                    sesh.tt_fname = t_files(jj).name;
                    
                    %% load spikes
                    disp('loading spikes')
                    spk_t = LoadSpikes({sesh.tt_fname});
                    spk_t = Data(spk_t{1});
                    
                    %% peri-event histograms
                    if PETH_generation == 1
                        sesh.PETH.Trial = 1; %generate PETH separated by cue conditions at trial start
                        sesh.PETH.Arm = 1; %generate PETH separated by arm location
                        sesh.PETH.Approach = 0; %generate PETH separated by approach and skip for unrewarded trials
                        sesh.PETH.Trial_np = 1; %generate PETH separated by cue conditions at trial start for when the rat approached
                        
                        PETH = genPETH(sesh,metadata,spk_t);
                    end
                    
                    %% rasterplots
                    RAST = genRAST(metadata,spk_t);
                    
                    %% extract event related firing rates
                    FRATE = genFRATE(metadata,spk_t);
                    
                    %% Wilcoxon signed rank test for cue-modulation
                    TESTS.WSR.Task.Trial_b4_vs_Trial = signrank(FRATE.Task.Trial_B4_firing_rate,FRATE.Task.Trial_firing_rate);
                    
                    %% save variables in .mat file
                    switch PETH_generation
                        case 0
                            save(cat(2,destination,sesh.session_id,'-TT',num2str(sesh.tt_number),'-cell',num2str(sesh.cell_number)),'metadata','FRATE','TESTS','RAST','spk_t','sesh');
                        case 1
                            save(cat(2,destination,sesh.session_id,'-TT',num2str(sesh.tt_number),'-cell',num2str(sesh.cell_number)),'metadata','FRATE','TESTS','RAST','PETH','spk_t','sesh');
                    end
                    clearvars -except f ii kj kk jj sesh rat_id day block_order_list directory destination metadata cell_count t_files str_length PETH_generation
                    cell_count = cell_count + 1;
                end
            end
        end
    end
end
end