% INPUTS:
% Data -- Standardized data you want to model. Rows are observations, 
% columns are features.
% y -- the dependent variable being predicted
% ShowPlot -- enter "1" if you want to see plots for the optimization, 
% enter "0" otherwise

% OUTPUTS:
% Mdl -- the optimized SVM model
% r -- the predicted values of the input y 

function [Mdl,r,imp] = tryEnsemble(Data,y,ShowPlot)
    
    hyperopts = ...
        struct('AcquisitionFunctionName','expected-improvement-plus',...
        'KFold',5,'ShowPlots',1,'Verbose',1);
    Mdl = fitrensemble(Data,y,...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',hyperopts);
    
    r = kfoldPredict(crossval(Mdl,'kfold',5));
    %r = predict(Mdl,Data);
    mn = floor(0.9*min([y';r])); mx = ceil(1.1*max([y';r]));
    %imp = [];
    imp = plotTreeImportance(Mdl,Data,y,ShowPlot);
    if ShowPlot == 1
        % plot the prediction
        figure; hold on;
        scatter(y,r);
        % plot(0:100,0:100,'r');
        plot(mn:mx,mn:mx,'r');
        xlabel('age (years)'); ylabel('predicted age (years)');
        title('Ensemble');

        rmse = sqrt(mean((y' - r).^2));
        % text(10,90,['rmse = ' num2str(rmse)])
        text(median(y),mx,['rmse = ' num2str(rmse)])
        xlim([mn mx]); ylim([mn mx]);
        
        
        
    end
    %}
    
end