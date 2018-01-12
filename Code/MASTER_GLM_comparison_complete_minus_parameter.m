%%
mat_files = dir('*.mat');

for kk = 1:length(dir('*.mat'))
    load(mat_files(kk).name);
    mat_overview.fname{kk} = mat_files(kk).name;
    disp(cat(2,num2str(kk),'/',num2str(length(dir('*.mat')))));
    
    for ll = 1:length(ALL_matrix(1,:))
        switch ll
            case 1
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.modality{kk} = removeTerms(mdl{kk},'Modality');
                    Comparison_GLM.Rsquared.modality(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.modality{kk}.Rsquared.Adjusted;
                end
            case 2
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.location{kk} = removeTerms(mdl{kk},'Location');
                    Comparison_GLM.Rsquared.location(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.location{kk}.Rsquared.Adjusted;
                end
            case 3
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.outcome{kk} = removeTerms(mdl{kk},'Outcome');
                    Comparison_GLM.Rsquared.outcome(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.outcome{kk}.Rsquared.Adjusted;
                end
            case 4
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.approach{kk} = removeTerms(mdl{kk},'Approach');
                    Comparison_GLM.Rsquared.approach(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.approach{kk}.Rsquared.Adjusted;
                end
            case 5
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.latency{kk} = removeTerms(mdl{kk},'Latency');
                    Comparison_GLM.Rsquared.latency(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.latency{kk}.Rsquared.Adjusted;
                end
            case 6
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.trial{kk} = removeTerms(mdl{kk},'Trial');
                    Comparison_GLM.Rsquared.trial(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.trial{kk}.Rsquared.Adjusted;
                end
            case 7
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.previous{kk} = removeTerms(mdl{kk},'Previous');
                    Comparison_GLM.Rsquared.previous(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.previous{kk}.Rsquared.Adjusted;
                end
            case 8
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.modxloc{kk} = removeTerms(mdl{kk},'Modality:Location');
                    Comparison_GLM.Rsquared.modxloc(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.modxloc{kk}.Rsquared.Adjusted;
                end
            case 9
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.modxout{kk} = removeTerms(mdl{kk},'Outcome:Modality');
                    Comparison_GLM.Rsquared.modxout(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.modxout{kk}.Rsquared.Adjusted;
                end
            case 10
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.locxout{kk} = removeTerms(mdl{kk},'Outcome:Location');
                    Comparison_GLM.Rsquared.locxout(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.locxout{kk}.Rsquared.Adjusted;
                end
            case 11
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.outxapp{kk} = removeTerms(mdl{kk},'Outcome:Approach');
                    Comparison_GLM.Rsquared.outxapp(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.outxapp{kk}.Rsquared.Adjusted;
                end
            case 12
                if ALL_matrix(kk,ll) == 1
                    Comparison_GLM.GLM.modxlocxout{kk} = removeTerms(mdl{kk},'Outcome:Modality:Location');
                    Comparison_GLM.Rsquared.modxlocxout(kk) = mdl{kk}.Rsquared.Adjusted - Comparison_GLM.GLM.modxlocxout{kk}.Rsquared.Adjusted;
                end
        end
    end
end

%%