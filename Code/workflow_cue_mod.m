mat_files = dir('*.mat');
cue_mod.all = 0;
cue_mod.FSI = 0;
cue_mod.MSN = 0;
cue_mod.Inc = 0;
cue_mod.Dec = 0;
cue_mod.MSN_Inc = 0;
cue_mod.MSN_Dec = 0;
cue_mod.FSI_Inc = 0;
cue_mod.FSI_Dec = 0;

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
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
    
    if RANK.two.Trial > 975 || RANK.two.Trial < 26
        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
            cue_mod.all = cue_mod.all + 1;
            
            switch class.all.sorted_isi{kk}(5) < 2
                                case 1
                                    cue_mod.FSI = cue_mod.FSI + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            cue_mod.FSI_Inc = cue_mod.FSI_Inc + 1;
                                            cue_mod.Inc = cue_mod.Inc + 1;
                                        case 0
                                            cue_mod.FSI_Dec = cue_mod.FSI_Dec + 1;
                                            cue_mod.Dec = cue_mod.Dec + 1;
                                    end
                                case 0
                                    cue_mod.MSN = cue_mod.MSN + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            cue_mod.MSN_Inc = cue_mod.MSN_Inc + 1;
                                            cue_mod.Inc = cue_mod.Inc + 1;
                                        case 0
                                            cue_mod.MSN_Dec = cue_mod.MSN_Dec + 1;
                                            cue_mod.Dec = cue_mod.Dec + 1;
                                    end
                            end
            
                    
            
        end
    end
            end
end