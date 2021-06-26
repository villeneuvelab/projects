clear
clc

yvar = 'age'; %leave this alone
testStr = []; %'DIAN'; % change if you want to look at a specific study (i.e. ADNI)
mut = []; %change this if you only want to look at mutation=1 or 0

% load chosen model and training results
load nn_search5.mat 

% load test data
load TestInfo_extra_new.mat

% get the data of interest
file = '/home/atbaria/Julie/Text_Files/data_combined_update_pibindex.xlsx';
%file = '/home/atbaria/Julie/Text_Files/PreventAD_extra.xlsx';
id = TestInfo.FeatureMatInfo.id;
if isempty(testStr)
    study = TestInfo.FeatureMatInfo.study; 
else
    study = testStr;
end
mutation = TestInfo.FeatureMatInfo.mutation;
[yTest,idxtest,idfin] = filter_y_data(file,id,study,yvar);
if ~isempty(mut)
    mutation = mutation(idxtest);
    idxtest = idxtest(mutation==mut); 
    yTest = yTest(mutation==mut);
end

% indicate the more interesting features
keepmetric = filterFeatures(TestInfo.names);

% get the graph metrics for the data chosen
DataTestRaw = real(TestInfo.GraphMetrics(idxtest,keepmetric));

% uncomment all these lines to remove outliers in the test set
outidxtest = removeOutliers(DataTestRaw,1,5);
DataTestRaw(outidxtest,:) = [];
yTest(outidxtest) = [];


% standardize test data by the train data mean and standard deviatio
[DataTest,~,~] = standardizeData(DataTestRaw,mean_train,sd_train);

% re-sort the features by importance
DataTestRank = DataTest(:,idxavg);

% test the network
ynntest = net(DataTestRank(:,1:numfeatures(numfeatureIdx))'); % model prediction
ypred = (conversionInfo.p(1)*yTest)+conversionInfo.p(2); % model fit
res = ynntest' - ypred; % the residual is the difference between the prediction and fit

% plot the model prediction only, against the actual age
plotScatterRmse_breakout(yTest,ynntest,...
    'predicted age',conversionInfo.p);
%xlim([10 100]); ylim([10 100])

% plot the model prediction, against the training fit
plotScatterMaturityIndex_breakout(yTest,...
    ynntest,miCnst,'predicted maturity',conversionInfo.pMI);
%xlim([10 100]); ylim([0 2.2])

% transform the predicted age by adjusting the slope of the training fit
% (residuals remain the same because transformation is linear)
x=0:120;
fx = (conversionInfo.p(1)*x) + conversionInfo.p(2);
pAdjust = polyfit(x,x-fx,1);
ynnAdjust = ynntest + (pAdjust(1)*yTest') + pAdjust(2);
figure; plot(x,x,'b');
hold on; plot(yTest,ynnAdjust,'r.');
xlabel('age'); ylabel('predicted age (adjusted)')
legend({'unity';'test data'});
xlim([10 100]); ylim([10 100])