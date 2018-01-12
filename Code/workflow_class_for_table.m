%%
mat_files = dir('*.mat');

summary_class.R053.HFN_Inc = 0;
summary_class.R053.HFN_Dec = 0;
summary_class.R053.SPN_Inc = 0;
summary_class.R053.SPN_Dec = 0;
summary_class.R056.HFN_Inc = 0;
summary_class.R056.HFN_Dec = 0;
summary_class.R056.SPN_Inc = 0;
summary_class.R056.SPN_Dec = 0;
summary_class.R057.HFN_Inc = 0;
summary_class.R057.HFN_Dec = 0;
summary_class.R057.SPN_Inc = 0;
summary_class.R057.SPN_Dec = 0;
summary_class.R060.HFN_Inc = 0;
summary_class.R060.HFN_Dec = 0;
summary_class.R060.SPN_Inc = 0;
summary_class.R060.SPN_Dec = 0;

summary_class.Cue.All = 0;
summary_class.Cue.HFN_Inc = 0;
summary_class.Cue.HFN_Dec = 0;
summary_class.Cue.SPN_Inc = 0;
summary_class.Cue.SPN_Dec = 0;

for kk = 1:length(dir('*.mat'))
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
    
    block_drift.block1_length(kk) = length(FRATE.Cue.Trial_firing_rate_block1);
    block_drift.block1_half(kk) = round(block_drift.block1_length(kk) / 2);
    block_drift.b1_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)));
    block_drift.b1_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    block_drift.MWU_b1(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block1(1:block_drift.block1_half(kk)),FRATE.Cue.Trial_firing_rate_block1(block_drift.block1_half(kk)+1:end));
    
    block_drift.block2_length(kk) = length(FRATE.Cue.Trial_firing_rate_block2);
    block_drift.block2_half(kk) = round(block_drift.block2_length(kk) / 2);
    block_drift.b2_1st_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)));
    block_drift.b2_2nd_avg(kk) = mean(FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    block_drift.MWU_b2(kk) = ranksum(FRATE.Cue.Trial_firing_rate_block2(1:block_drift.block2_half(kk)),FRATE.Cue.Trial_firing_rate_block2(block_drift.block2_half(kk)+1:end));
    
    class.all.isi{kk} = diff(spk_t);
    class.all.median_isi(kk) = median(class.all.isi{kk});
    class.all.frate(kk) = FRATE.Overall.firing_rate_total;
    class.all.sorted_isi{kk} = sort(class.all.isi{kk},'descend');
    
    switch sesh.session_id(1:4)
        case 'R053'
            switch class.all.sorted_isi{kk}(5) < 2
                case 1
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R053.HFN_Inc = summary_class.R053.HFN_Inc + 1;
                        case 0
                            summary_class.R053.HFN_Dec = summary_class.R053.HFN_Dec + 1;
                    end
                case 0
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R053.SPN_Inc = summary_class.R053.SPN_Inc + 1;
                        case 0
                            summary_class.R053.SPN_Dec = summary_class.R053.SPN_Dec + 1;
                    end
            end
            case 'R056'
            switch class.all.sorted_isi{kk}(5) < 2
                case 1
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R056.HFN_Inc = summary_class.R056.HFN_Inc + 1;
                        case 0
                            summary_class.R056.HFN_Dec = summary_class.R056.HFN_Dec + 1;
                    end
                case 0
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R056.SPN_Inc = summary_class.R056.SPN_Inc + 1;
                        case 0
                            summary_class.R056.SPN_Dec = summary_class.R056.SPN_Dec + 1;
                    end
            end
            case 'R057'
            switch class.all.sorted_isi{kk}(5) < 2
                case 1
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R057.HFN_Inc = summary_class.R057.HFN_Inc + 1;
                        case 0
                            summary_class.R057.HFN_Dec = summary_class.R057.HFN_Dec + 1;
                    end
                case 0
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R057.SPN_Inc = summary_class.R057.SPN_Inc + 1;
                        case 0
                            summary_class.R057.SPN_Dec = summary_class.R057.SPN_Dec + 1;
                    end
            end
            case 'R060'
            switch class.all.sorted_isi{kk}(5) < 2
                case 1
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R060.HFN_Inc = summary_class.R060.HFN_Inc + 1;
                        case 0
                            summary_class.R060.HFN_Dec = summary_class.R060.HFN_Dec + 1;
                    end
                case 0
                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                        case 1
                            summary_class.R060.SPN_Inc = summary_class.R060.SPN_Inc + 1;
                        case 0
                            summary_class.R060.SPN_Dec = summary_class.R060.SPN_Dec + 1;
                    end
            end
    end
    
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
                    if RANK.two.Trial > 975 || RANK.two.Trial < 26
                        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                            summary_class.Cue.All = summary_class.Cue.All + 1;
                            switch class.all.sorted_isi{kk}(5) < 2
                                case 1
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            summary_class.Cue.HFN_Inc = summary_class.Cue.HFN_Inc + 1;
                                        case 0
                                            summary_class.Cue.HFN_Dec = summary_class.Cue.HFN_Dec + 1;
                                    end
                                case 0
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            summary_class.Cue.SPN_Inc = summary_class.Cue.SPN_Inc + 1;
                                        case 0
                                            summary_class.Cue.SPN_Dec = summary_class.Cue.SPN_Dec + 1;
                                    end
                            end
                        end
                    end
            end
end