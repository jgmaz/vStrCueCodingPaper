%% cue onset vs NP

for iFeature = 1:3
    for iCell = 1:length(ALL_matrix)
        if ALL_matrix(iCell,iFeature) == 1 && ALL_matrix_NP(iCell,iFeature) == 1
            overlap.cells(iCell,iFeature) = 1;
        else
            overlap.cells(iCell,iFeature) = 0;
        end
    end
    
    overlap.total(1,iFeature) = sum(overlap.cells(:,iFeature));
    overlap.percentage(1,iFeature) = overlap.total(1,iFeature)/sum(ALL_matrix_NP(:,iFeature));
    prop1 = sum(ALL_matrix(:,iFeature)) / 133;
        prop2 = sum(ALL_matrix_NP(:,iFeature)) / 133;
        overlap.exp(1,iFeature) = prop1*prop2*133;
        overlap.chi(1,iFeature) = (overlap.exp(1,iFeature) - overlap.total(1,iFeature))^2 / overlap.exp(1,iFeature);
        overlap.pvalue(1,iFeature) = 1 - chi2cdf(overlap.chi(1,iFeature),1);
       
end

%% for cue-onset GLM
Variables = {'Mod' 'Loc' 'Out' 'App' 'Lat' 'Trial' 'Prev'};
for iVar1 = 1:3%length(Variables)
    for iVar2 = 1:3%length(Variables)
        Overlap.obs.(Variables{iVar1}).(Variables{iVar2}) = 0;
        for iCell = 1:length(ALL_matrix)
           if ALL_matrix(iCell,iVar1) == 1 && ALL_matrix(iCell,iVar2) == 1
        Overlap.obs.(Variables{iVar1}).(Variables{iVar2}) = Overlap.obs.(Variables{iVar1}).(Variables{iVar2}) + 1;
           end 
        end
        prop1 = sum(ALL_matrix(:,iVar1)) / 133;
        prop2 = sum(ALL_matrix(:,iVar2)) / 133;
        Overlap.exp.(Variables{iVar1}).(Variables{iVar2}) = prop1*prop2*133;
        Overlap.chi.(Variables{iVar1}).(Variables{iVar2}) = (Overlap.obs.(Variables{iVar1}).(Variables{iVar2}) -...
            Overlap.exp.(Variables{iVar1}).(Variables{iVar2}))^2 / Overlap.exp.(Variables{iVar1}).(Variables{iVar2});
        Overlap.pvalue.(Variables{iVar1}).(Variables{iVar2}) = 1 - chi2cdf(Overlap.chi.(Variables{iVar1}).(Variables{iVar2}),1);
        if Overlap.pvalue.(Variables{iVar1}).(Variables{iVar2}) < .05
            disp(cat(2,Variables(iVar1),Variables(iVar2)));
        end
    end
end
%%
figure
venn([sum(ALL_matrix(:,1)) sum(ALL_matrix(:,2)) sum(ALL_matrix(:,3))],...
    [Overlap.ModxLoc Overlap.ModxOut Overlap.LocxOut Overlap.ModxLocxOut],'ErrMinMode','TotalError');
