%%
mat_files = dir('*.mat');

Predictors = {'Modality' 'Location' 'Outcome'};
count_cue = 1;

for iPred = 1:length(Predictors)
     GLM_coeff.(Predictors{iPred}).ALL = zeros(1,133);
    for iPred2 = 1:length(Predictors)
        if iPred ~= iPred2
            count.(Predictors{iPred}).(Predictors{iPred2}) = 1;
        end
    end
end

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
    
    switch isempty(mdl{kk})
        case 0
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
                    if RANK.two.Trial > 975 || RANK.two.Trial < 26
                        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                            for ll = 2:length(mdl{1,kk}.CoefficientNames)
                                if mdl{1,kk}.Coefficients{ll,4} < .01 %&& mdl{1,kk}.Coefficients{ll,4} ~= 0
                                    for iPred = 1:length(Predictors)
                                        if strcmp(mdl{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                            GLM_coeff.(Predictors{iPred}).ALL(count_cue) = 1;
                                            for nn = 2:length(mdl{1,kk}.CoefficientNames)
                                                if mdl{1,kk}.Coefficients{nn,4} < .01
                                                for iPred2 = 1:length(Predictors)
                                                    if iPred ~= iPred2
                                                        if strcmp(mdl{1,kk}.CoefficientNames(nn),Predictors(iPred2))
                                                            GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred})(count.(Predictors{iPred}).(Predictors{iPred2})) = mdl{1,kk}.Coefficients{ll,1};
                                                            GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2})(count.(Predictors{iPred}).(Predictors{iPred2})) = mdl{1,kk}.Coefficients{nn,1};
                                                            count.(Predictors{iPred}).(Predictors{iPred2}) = count.(Predictors{iPred}).(Predictors{iPred2}) + 1;
                                                        end
                                                    end
                                                end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            count_cue = count_cue + 1;
                        end
                    end
            end
    end
end

%% Coerr coeff

for iPred = 1:length(Predictors)
    for iPred2 = 1:length(Predictors)
        if iPred ~= iPred2
            [GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).AbsCorr GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).AbsCorrPvalues] = corrcoef(abs(GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred})),abs(GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2})))
            [GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).RawCorr GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).RawCorrPvalues] = corrcoef(GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred}),GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2}))
            [GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr GLM_coeff.(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues] = corrcoef( GLM_coeff.(Predictors{iPred}).ALL, GLM_coeff.(Predictors{iPred2}).ALL)
        end
    end
end

%% across epochs
mat_files = dir('*.mat');

Predictors = {'Modality' 'Location' 'Outcome'};
count_cue = 1;

for iPred = 1:length(Predictors)
     GLM_coeff.(Predictors{iPred}).ALL.epoch1 = zeros(1,133);
     GLM_coeff.(Predictors{iPred}).ALL.epoch2 = zeros(1,133);
            count.(Predictors{iPred}) = 1;

end

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
    
    switch isempty(mdl{kk})
        case 0
            switch isempty(mdl_epoch2{kk})
                case 0
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
                    if RANK.two.Trial > 975 || RANK.two.Trial < 26
                        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                            for ll = 2:length(mdl{1,kk}.CoefficientNames)
                                if mdl{1,kk}.Coefficients{ll,4} < .01 %&& mdl{1,kk}.Coefficients{ll,4} ~= 0
                                    for iPred = 1:length(Predictors)
                                        if strcmp(mdl{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                            GLM_coeff.(Predictors{iPred}).ALL.epoch1(count_cue) = 1;
                                            for nn = 2:length(mdl_epoch2{1,kk}.CoefficientNames)
                                                if mdl_epoch2{1,kk}.Coefficients{nn,4} < .01
                                                        if strcmp(mdl_epoch2{1,kk}.CoefficientNames(nn),Predictors(iPred))
                                                            GLM_coeff.(Predictors{iPred}).epoch1(count.(Predictors{iPred})) = mdl{1,kk}.Coefficients{ll,1};
                                                            GLM_coeff.(Predictors{iPred}).epoch2(count.(Predictors{iPred})) = mdl_epoch2{1,kk}.Coefficients{nn,1};
                                                            count.(Predictors{iPred}) = count.(Predictors{iPred}) + 1;
                                                        end

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            for ll = 2:length(mdl_epoch2{1,kk}.CoefficientNames)
                                if mdl_epoch2{1,kk}.Coefficients{ll,4} < .01 %&& mdl{1,kk}.Coefficients{ll,4} ~= 0
                                    for iPred = 1:length(Predictors)
                                        if strcmp(mdl_epoch2{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                            GLM_coeff.(Predictors{iPred}).ALL.epoch2(count_cue) = 1;
                                        end
                                    end
                                end
                            end
                            count_cue = count_cue + 1;
                        end
                    end
            end
            end
    end
end

%% Coerr coeff

for iPred = 1:length(Predictors)
            [GLM_coeff.(Predictors{iPred}).AbsCorr GLM_coeff.(Predictors{iPred}).AbsCorrPvalues] = corrcoef(abs(GLM_coeff.(Predictors{iPred}).epoch1),abs(GLM_coeff.(Predictors{iPred}).epoch2))
            [GLM_coeff.(Predictors{iPred}).RawCorr GLM_coeff.(Predictors{iPred}).RawCorrPvalues] = corrcoef(GLM_coeff.(Predictors{iPred}).epoch1,GLM_coeff.(Predictors{iPred}).epoch2)
            [GLM_coeff.(Predictors{iPred}).RecodeCorr GLM_coeff.(Predictors{iPred}).RecodeCorrPvalues] = corrcoef(GLM_coeff.(Predictors{iPred}).ALL.epoch1,GLM_coeff.(Predictors{iPred}).ALL.epoch2)
end