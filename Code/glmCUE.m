function yfit = glmCUE(xtrain,ytrain,xtest)

model = fitglm(xtrain,ytrain,'Link','logit','Distribution','binomial');

switch length(xtest(1,:))
    case 1
        for iTrial = 1:length(xtest)
            yfit(iTrial) = exp(model.Coefficients{1,1} + (xtest(iTrial) * model.Coefficients{2,1})) / (1 + exp(model.Coefficients{1,1} + (xtest(iTrial) * model.Coefficients{2,1})));
        end
    case 2
        for iTrial = 1:length(xtest)
            yfit(iTrial) = exp(model.Coefficients{1,1} + (xtest(iTrial,1) * model.Coefficients{2,1}) + (xtest(iTrial,2) * model.Coefficients{3,1})) / (1 + exp(model.Coefficients{1,1} + (xtest(iTrial,1) * model.Coefficients{2,1}) + (xtest(iTrial,2) * model.Coefficients{3,1})));
        end
end

yfit(isnan(yfit)) = 0;