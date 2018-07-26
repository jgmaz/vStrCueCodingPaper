%%
mat_files = dir('*.mat');

Predictors = {'Modality' 'Location' 'Outcome'};
count_cue = 1;

for iPred = 1:length(Predictors)
    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL = zeros(1,133);
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
    
    switch isempty(mdl_epoch1{kk})
        case 0
            switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                case 0
                    if RANK.two.Trial > 975 || RANK.two.Trial < 26
                        if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                            for ll = 2:length(mdl_epoch1{1,kk}.CoefficientNames)
                                if mdl_epoch1{1,kk}.Coefficients{ll,4} < .01 %&& mdl_epoch1{1,kk}.Coefficients{ll,4} ~= 0
                                    for iPred = 1:length(Predictors)
                                        if strcmp(mdl_epoch1{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                            GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL(count_cue) = 1;
                                            for nn = 2:length(mdl_epoch1{1,kk}.CoefficientNames)
                                                if mdl_epoch1{1,kk}.Coefficients{nn,4} < .01
                                                    for iPred2 = 1:length(Predictors)
                                                        if iPred ~= iPred2
                                                            if strcmp(mdl_epoch1{1,kk}.CoefficientNames(nn),Predictors(iPred2))
                                                                GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred})(count.(Predictors{iPred}).(Predictors{iPred2})) = mdl_epoch1{1,kk}.Coefficients{ll,1};
                                                                GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2})(count.(Predictors{iPred}).(Predictors{iPred2})) = mdl_epoch1{1,kk}.Coefficients{nn,1};
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
            [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).AbsCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).AbsCorrPvalues] = corrcoef(abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred})),abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2})))
            [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).RawCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).RawCorrPvalues] = corrcoef(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred}),GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).(Predictors{iPred2}))
            [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).(Predictors{iPred2}).RecodeCorrPvalues] = corrcoef( GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL, GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred2}).ALL)
        end
    end
end

%% across epochs
Epoch = {'cueon' 'NP' 'outcome' 'cueoff'}; %1 = cue on, 2 = NP, 3 = outcome, 4 = cue off
for iEpoch = 1:3%:length(Epoch)
    disp(iEpoch)
    switch iEpoch
        case 1
            % mat_files = dir('2018-03-12-GLM_cueon*');
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous' ...
            %     'ModxLoc' 'ModxOut' 'LocxOut' 'OutxApp' 'ModxLocxOut'};
            % file_order = {'2018-05-24-GLM_cueon_-0.5-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.4-recodePrevTrial.mat' '2018-05-24_cueon_-0.3-recodePrevTrial.mat' ...
            %     '2018-05-24-GLM_cueon_-0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.1-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.1-recodePrevTrial.mat' ...
            %     '2018-05-24-GLM_cueon_0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.3-recodePrevTrial.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
            file_order = {'2018-03-12-GLM_cueon_-0.5.mat' '2018-03-12-GLM_cueon_-0.4.mat' '2018-03-12-GLM_cueon_-0.3.mat' ...
                '2018-03-12-GLM_cueon_-0.2.mat' '2018-03-12-GLM_cueon_-0.1.mat' '2018-03-12-GLM_cueon_0.mat' '2018-03-12-GLM_cueon_0.1.mat' ...
                '2018-03-12-GLM_cueon_0.2.mat' '2018-03-12-GLM_cueon_0.3.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
        case 2
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-12-GLM_NP_-0.5-500bin.mat' '2018-03-12-GLM_NP_-0.4-500bin.mat' '2018-03-12-GLM_NP_-0.3-500bin.mat' ...
                '2018-03-12-GLM_NP_-0.2-500bin.mat' '2018-03-12-GLM_NP_-0.1-500bin.mat' '2018-03-12-GLM_NP_0-500bin.mat' '2018-03-12-GLM_NP_0.1-500bin.mat' ...
                '2018-03-12-GLM_NP_0.2-500bin.mat' '2018-03-12-GLM_NP_0.3-500bin.mat' '2018-03-12-GLM_NP_0.4-500bin.mat' '2018-03-12-GLM_NP_0.5-500bin.mat'};
        case 3
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-23-GLM_outcome_-0.5.mat' '2018-03-23-GLM_outcome_-0.4.mat' '2018-03-23-GLM_outcome_-0.3.mat' ...
                '2018-03-23-GLM_outcome_-0.2.mat' '2018-03-23-GLM_outcome_-0.1.mat' '2018-03-23-GLM_outcome_0.mat' '2018-03-23-GLM_outcome_0.1.mat' ...
                '2018-03-23-GLM_outcome_0.2.mat' '2018-03-23-GLM_outcome_0.3.mat' '2018-03-23-GLM_outcome_0.4.mat' '2018-03-23-GLM_outcome_0.5.mat'};
        case 4
            mat_files = dir('E:\Jimmie\Jimmie\Analysis\2018-03-12-GLM_cueoff*');
            % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
            file_order = {'2018-03-12-GLM_cueoff_-0.5.mat' '2018-03-12-GLM_cueoff_-0.4.mat' '2018-03-12-GLM_cueoff_-0.3.mat' ...
                '2018-03-12-GLM_cueoff_-0.2.mat' '2018-03-12-GLM_cueoff_-0.1.mat' '2018-03-12-GLM_cueoff_0.mat' '2018-03-12-GLM_cueoff_0.1.mat' ...
                '2018-03-12-GLM_cueoff_0.2.mat' '2018-03-12-GLM_cueoff_0.3.mat' '2018-03-12-GLM_cueoff_0.4.mat' '2018-03-12-GLM_cueoff_0.5.mat'};
    end
    
    for iEpoch2 = 1:3
        disp(iEpoch2)
        switch iEpoch2
            case 1
                % mat_files = dir('2018-03-12-GLM_cueon*');
                mat_files2 = dir('E:\Jimmie\Jimmie\Analysis\2018-03*');
                % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome' 'Approach' 'Latency' 'Trial' 'Previous' ...
                %     'ModxLoc' 'ModxOut' 'LocxOut' 'OutxApp' 'ModxLocxOut'};
                % file_order = {'2018-05-24-GLM_cueon_-0.5-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.4-recodePrevTrial.mat' '2018-05-24_cueon_-0.3-recodePrevTrial.mat' ...
                %     '2018-05-24-GLM_cueon_-0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_-0.1-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.1-recodePrevTrial.mat' ...
                %     '2018-05-24-GLM_cueon_0.2-recodePrevTrial.mat' '2018-05-24-GLM_cueon_0.3-recodePrevTrial.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
                file_order2 = {'2018-03-12-GLM_cueon_-0.5.mat' '2018-03-12-GLM_cueon_-0.4.mat' '2018-03-12-GLM_cueon_-0.3.mat' ...
                    '2018-03-12-GLM_cueon_-0.2.mat' '2018-03-12-GLM_cueon_-0.1.mat' '2018-03-12-GLM_cueon_0.mat' '2018-03-12-GLM_cueon_0.1.mat' ...
                    '2018-03-12-GLM_cueon_0.2.mat' '2018-03-12-GLM_cueon_0.3.mat' '2018-03-12-GLM_cueon_0.4.mat' '2018-03-12-GLM_cueon_0.5.mat'};
            case 2
                mat_files2 = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
                % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
                file_order2 = {'2018-03-12-GLM_NP_-0.5-500bin.mat' '2018-03-12-GLM_NP_-0.4-500bin.mat' '2018-03-12-GLM_NP_-0.3-500bin.mat' ...
                    '2018-03-12-GLM_NP_-0.2-500bin.mat' '2018-03-12-GLM_NP_-0.1-500bin.mat' '2018-03-12-GLM_NP_0-500bin.mat' '2018-03-12-GLM_NP_0.1-500bin.mat' ...
                    '2018-03-12-GLM_NP_0.2-500bin.mat' '2018-03-12-GLM_NP_0.3-500bin.mat' '2018-03-12-GLM_NP_0.4-500bin.mat' '2018-03-12-GLM_NP_0.5-500bin.mat'};
            case 3
                mat_files2 = dir('E:\Jimmie\Jimmie\Analysis\2018-03-*');
                % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
                file_order2 = {'2018-03-23-GLM_outcome_-0.5.mat' '2018-03-23-GLM_outcome_-0.4.mat' '2018-03-23-GLM_outcome_-0.3.mat' ...
                    '2018-03-23-GLM_outcome_-0.2.mat' '2018-03-23-GLM_outcome_-0.1.mat' '2018-03-23-GLM_outcome_0.mat' '2018-03-23-GLM_outcome_0.1.mat' ...
                    '2018-03-23-GLM_outcome_0.2.mat' '2018-03-23-GLM_outcome_0.3.mat' '2018-03-23-GLM_outcome_0.4.mat' '2018-03-23-GLM_outcome_0.5.mat'};
            case 4
                mat_files2 = dir('E:\Jimmie\Jimmie\Analysis\2018-03-12-GLM_cueoff*');
                % predictor_list = {'Cue' 'Modality' 'Location' 'Outcome'};
                file_order2 = {'2018-03-12-GLM_cueoff_-0.5.mat' '2018-03-12-GLM_cueoff_-0.4.mat' '2018-03-12-GLM_cueoff_-0.3.mat' ...
                    '2018-03-12-GLM_cueoff_-0.2.mat' '2018-03-12-GLM_cueoff_-0.1.mat' '2018-03-12-GLM_cueoff_0.mat' '2018-03-12-GLM_cueoff_0.1.mat' ...
                    '2018-03-12-GLM_cueoff_0.2.mat' '2018-03-12-GLM_cueoff_0.3.mat' '2018-03-12-GLM_cueoff_0.4.mat' '2018-03-12-GLM_cueoff_0.5.mat'};
        end
        for iGLM = 1:length(file_order)
            for iFind = 1:length(mat_files)
                if strcmp(mat_files(iFind).name,file_order{iGLM}) == 1
                    load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files(iFind).name),'mdl');
                    mdl_epoch1 = mdl;
                    
                    for iGLM2 = 1:length(file_order2)
                        for iFind2 = 1:length(mat_files2)
                            if strcmp(mat_files2(iFind2).name,file_order2{iGLM2}) == 1
                                load(strcat('E:\Jimmie\Jimmie\Analysis\',mat_files2(iFind2).name),'mdl');
                                mdl_epoch2 = mdl;
                                %                 disp(cat(2,num2str(iGLM),'/',num2str(length(mat_files))));
                                
                                
                                cell_files = dir('*.mat');
                                
                                Predictors = {'Modality' 'Location' 'Outcome'};
                                count_cue = 1;
                                
                                for iPred = 1:length(Predictors)
                                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch1 = zeros(1,133);
                                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch2 = zeros(1,133);
                                    count.(Predictors{iPred}) = 1;
                                    
                                end
                                
                                for kk = 1:length(dir('*.mat'))
                                    load(cell_files(kk).name);
                                    mat_overview.fname{kk} = cell_files(kk).name;
%                                     disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat'))),' for epochs ',num2str(iEpoch),' (file #',num2str(iGLM),') & ',num2str(iEpoch2),' (file #',num2str(iGLM2),')'));
                                    
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
                                    
                                    switch isempty(mdl_epoch1{kk})
                                        case 0
                                            switch isempty(mdl_epoch2{kk})
                                                case 0
                                                    switch block_drift.MWU_b1(kk) < .01 || block_drift.MWU_b2(kk) < .01
                                                        case 0
                                                            if RANK.two.Trial > 975 || RANK.two.Trial < 26
                                                                if TESTS.WSR.Task.Trial_b4_vs_Trial < .01
                                                                    disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat'))),' for epochs ',num2str(iEpoch),' (file #',num2str(iGLM),') & ',num2str(iEpoch2),' (file #',num2str(iGLM2),')'));
                                                                    for ll = 2:length(mdl_epoch1{1,kk}.CoefficientNames)
                                                                        if mdl_epoch1{1,kk}.Coefficients{ll,4} < .01 %&& mdl_epoch1{1,kk}.Coefficients{ll,4} ~= 0
                                                                            for iPred = 1:length(Predictors)
                                                                                if strcmp(mdl_epoch1{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                                                                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch1(count_cue) = 1;
                                                                                    for nn = 2:length(mdl_epoch2{1,kk}.CoefficientNames)
                                                                                        if mdl_epoch2{1,kk}.Coefficients{nn,4} < .01
                                                                                            if strcmp(mdl_epoch2{1,kk}.CoefficientNames(nn),Predictors(iPred))
                                                                                                GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch1(count.(Predictors{iPred})) = mdl_epoch1{1,kk}.Coefficients{ll,1};
                                                                                                GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch2(count.(Predictors{iPred})) = mdl_epoch2{1,kk}.Coefficients{nn,1};
                                                                                                count.(Predictors{iPred}) = count.(Predictors{iPred}) + 1;
                                                                                            end
                                                                                            
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                    for ll = 2:length(mdl_epoch2{1,kk}.CoefficientNames)
                                                                        if mdl_epoch2{1,kk}.Coefficients{ll,4} < .01 %&& mdl_epoch1{1,kk}.Coefficients{ll,4} ~= 0
                                                                            for iPred = 1:length(Predictors)
                                                                                if strcmp(mdl_epoch2{1,kk}.CoefficientNames(ll),Predictors(iPred))
                                                                                    GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch2(count_cue) = 1;
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
                                for iPred = 1:length(Predictors)
%                                     [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).AbsCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).AbsCorrPvalues] = corrcoef(abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch1),abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch2));
%                                     [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RawCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RawCorrPvalues] = corrcoef(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch1,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch2);
                                    [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues] = corrcoef(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch1,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch2);
                                    
                                    GLM_coeffs.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(iGLM,iGLM2) = GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr(2);
                                    GLM_coeffs.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(iGLM,iGLM2) = GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues(2);
                                end
                                GLM_coeff = [];
                            end
                        end
                    end
                end
            end
        end
    end
end
%% Coerr coeff

% for iPred = 1:length(Predictors)
%     [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).AbsCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).AbsCorrPvalues] = corrcoef(abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch1),abs(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch2))
%     [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RawCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RawCorrPvalues] = corrcoef(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch1,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).epoch2)
%     [GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorr GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).RecodeCorrPvalues] = corrcoef(GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch1,GLM_coeff.(Epoch{iEpoch}).(Epoch{iEpoch2}).(Predictors{iPred}).ALL.epoch2)
% end