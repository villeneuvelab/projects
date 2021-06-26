

function [net,y] = tryNN(Data0,y0,hiddenLayerSize,showPlot)
    
    %setdemorandstream(491218382);

    % Choose a Training Function
    % For a list of all training functions type: help nntrain
    % 'trainlm' is usually fastest.
    % 'trainbr' takes longer but may be better for challenging problems.
    % 'trainscg' uses less memory. Suitable in low memory situations.
    trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.

    t = y0';
    x = Data0';
    
    net = fitnet(hiddenLayerSize,trainFcn);

    % Choose Input and Output Pre/Post-Processing Functions
    % For a list of all processing functions type: help nnprocess
    net.input.processFcns = {'removeconstantrows','mapminmax'};
    net.output.processFcns = {'removeconstantrows','mapminmax'};

    % Setup Division of Data for Training, Validation, Testing
    % For a list of all data division functions type: help nndivide
    net.divideFcn = 'dividerand';  % Divide data interleaved
    net.divideMode = 'sample';  % Divide up every sample
    net.divideParam.trainRatio = 80/100;
    net.divideParam.valRatio = 10/100;
    net.divideParam.testRatio = 10/100;
    
    % other options
    net.trainParam.showWindow = 0; % do not show gui

    % Choose a Performance Function
    % For a list of all performance functions type: help nnperformance
    net.performFcn = 'mse';  % Mean Squared Error

    % Choose Plot Functions
    % For a list of all plot functions type: help nnplot
    net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
        'plotregression', 'plotfit'};

    % Train the Network
    [net,tr] = train(net,x,t);

    % Test the Network
    y = net(x);
    e = gsubtract(t,y);
    %performance = perform(net,t,y);

%         % Recalculate Training, Validation and Test Performance
%         trainTargets = t .* tr.trainMask{1};
%         valTargets = t .* tr.valMask{1};
%         testTargets = t .* tr.testMask{1};
%         trainPerformance = perform(net,trainTargets,y);
%         valPerformance = perform(net,valTargets,y);
%         testPerformance = perform(net,testTargets,y);

        % View the Network
        %view(net)
        
    
    if showPlot == 1
        plotScatterRmse(t,y,'NeuralNetwork');
    end
        
        
        % Plots
        % Uncomment these lines to enable various plots.
        %figure, plotperform(tr)
        %figure, plottrainstate(tr)
        %figure, ploterrhist(e)
        %figure, plotregression(t,y)
        %figure, plotfit(net,x,t)
end

   