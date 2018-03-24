block_drift.HFN_Inc = 0;
block_drift.SPN_Inc = 0;
block_drift.Inc = 0;
 block_drift.HFN = 0;
  block_drift.SPN = 0;

block_drift.HFN_Dec = 0;
block_drift.SPN_Dec = 0;
block_drift.Dec = 0;

mat_files = dir('*.mat');

for kk = 1:length(dir('*.mat'))
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
    
switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
    
    class.all.isi{kk} = diff(spk_t);
    class.all.median_isi(kk) = median(class.all.isi{kk});
    class.all.frate(kk) = FRATE.Overall.firing_rate_total;
    class.all.sorted_isi{kk} = sort(class.all.isi{kk},'descend');
    
    switch class.all.sorted_isi{kk}(5) < 2
                                case 1
                                    block_drift.HFN = block_drift.HFN + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            block_drift.HFN_Inc = block_drift.HFN_Inc + 1;
                                            block_drift.Inc = block_drift.Inc + 1;
                                        case 0
                                            block_drift.HFN_Dec = block_drift.HFN_Dec + 1;
                                            block_drift.Dec = block_drift.Dec + 1;
                                    end
                                case 0
                                    block_drift.SPN = block_drift.SPN + 1;
                                    switch mean(FRATE.Task.Trial_firing_rate) > mean(FRATE.Task.Trial_B4_firing_rate)
                                        case 1
                                            block_drift.SPN_Inc = block_drift.SPN_Inc + 1;
                                            block_drift.Inc = block_drift.Inc + 1;
                                        case 0
                                            block_drift.SPN_Dec = block_drift.SPN_Dec + 1;
                                            block_drift.Dec = block_drift.Dec + 1;
                                    end
    end
                            
end
end