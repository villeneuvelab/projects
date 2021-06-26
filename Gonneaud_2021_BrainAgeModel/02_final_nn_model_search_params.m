
clear
clc

yvar = 'age'; % leave this alone
testStr = 'ICBM'; % ICBM was always used as the validation set  

% load train data
load TrainInfo_extra_new.mat

% filter data by the available dependent variable
file = '/home/atbaria/Julie/Text_Files/data_combined_update_pibindex.xlsx';
id = TrainInfo.FeatureMatInfo.id;
study = TrainInfo.FeatureMatInfo.study;
[y,finidx_train] = filter_y_data(file,id,study,yvar);
DataRaw = real(TrainInfo.GraphMetrics(finidx_train,:));

%variable names (graph metrics)
names = TrainInfo.names;

% remove perhaps unimportant metrics (determined from visual inspection)
keepmetric = []; ct = 0;
for i = 1:length(names)
    if isempty([strfind(names{i},'iqr'), ...
            strfind(names{i},'25'), strfind(names{i},'75')])
        ct = ct+1;
        keepmetric(ct) = i;
    end
end
DataRaw = DataRaw(:,keepmetric);

% remove outliers
outidxtrain = removeOutliers(DataRaw,1,5);
DataRaw(outidxtrain,:) = [];
y(outidxtrain) = [];

% standardize the data across examples
[Data,mean_train,sd_train] = standardizeData(DataRaw,[],[]);

% load the test data and standardize to the train set
load TestInfo_extra_new.mat

% filter by y-variable
id = TestInfo.FeatureMatInfo.id;
study = TestInfo.FeatureMatInfo.study;
[yTest,finidx_temp] = filter_y_data(file,id,study,yvar);

% if filtering further by study
if ~isempty(testStr)
    testStudyIdx = find(strcmp(TestInfo.FeatureMatInfo.study,testStr));
    [testStudyIdx,ia,~] = intersect(finidx_temp,testStudyIdx);
    yTest = yTest(ia);
else
    testStudyIdx = finidx_temp;
end

DataTestRaw = real(TestInfo.GraphMetrics(testStudyIdx,keepmetric));

% remove outliers in the test set
outidxtest = removeOutliers(DataTestRaw,1,5);
DataTestRaw(outidxtest,:) = [];
yTest(outidxtest) = [];

% standardize test data by the train data mean and standard deviatio
[DataTest,~,~] = standardizeData(DataTestRaw,mean_train,sd_train);

% build a couple models to determine feature importance
[MdlTree,rtree,impTree] = tryEnsemble(DataRaw,y',0);
[MdlSVM,rsvm] = trySVMLinear(Data,y',0);

% sort features by importance
impSVM = abs(MdlSVM.Beta);
impSVMNorm = (impSVM-min(impSVM))/(max(impSVM)-min(impSVM));
impTreeNorm = (impTree-min(impTree))/(max(impTree)-min(impTree));
figure; plot(impTreeNorm,impSVMNorm,'r.','MarkerSize',30); 
text(impTreeNorm,impSVMNorm,names(keepmetric)); xlabel('ensemble importance');
ylabel('svm importance'); xlim([-0.1 1.1]); ylim([-0.1 1.1]);

% perform ranking of features
n = length(impSVM);
rankv = 0:1/n:1; rankv = flipdim(rankv,2);
[~,idxsvm] = sort(impSVM,'descend'); 
ranksvm = zeros(1,n); for i = 1:n; ranksvm(idxsvm(i)) = rankv(i); end
[~,idxtree] = sort(impTree,'descend'); 
ranktree = zeros(1,n); for i = 1:n; ranktree(idxtree(i)) = rankv(i); end
rankavg = (ranktree + ranksvm)/2;
[~,idxavg] = sort(rankavg,'descend');
DataRank = Data(:,idxavg);
DataTestRank = DataTest(:,idxavg);

% Now test ranked features in multiple neural network models
[PERFORMANCE, NETS, numfeatures] = ...
    nnSearchParams(DataRank,DataTestRank,y,yTest);
%}
% plot the performance of the neural network models
plotNNSearch(PERFORMANCE.rmse_train,numfeatures,NETS(1,:),'training data');
plotNNSearch(PERFORMANCE.rmse_test,numfeatures,NETS(1,:),'ICBM data');
plotNNSearch(PERFORMANCE.rmse_null,numfeatures,NETS(1,:),'null');

% find the model with the best performance, and show the scatter plots of
% the predicted age against the actual age
rmin = PERFORMANCE.rmse_test;
[numfeatureIdx,c] = find(rmin==min(rmin(:)));
units = NETS{numfeatureIdx,c}; % number of features, network
[net,ynntrain] = tryNN(DataRank(:,1:numfeatures(numfeatureIdx)),y,units,0);
ynntest = net(DataTestRank(:,1:numfeatures(numfeatureIdx))');
plotScatterRmse(y,ynntrain,'training data')
plotScatterRmse(yTest,ynntest,'ICBM data')

% you can also get some index to address the bias here, like a maturity
% index, by judging how well the test data fits to the training model
% (rather than the actual age)
conversionInfo = []; 
p = polyfit(y,ynntrain',1); 
conversionInfo.p = p; % linear fit of age to predicted age, training

miCnst = mean(ynntrain);
pMI = polyfit(y,ynntrain'/miCnst,1);
conversionInfo.pMI = pMI; % linear fit of age to maturity index, training
figure; hold on; 
plotScatterMaturityIndex(y,ynntrain/miCnst,'train',conversionInfo.pMI,'k.',3);
plotScatterMaturityIndex(yTest,ynntest/miCnst,'ICBM',conversionInfo.pMI,'r.',30);

% SAVE THIS STUFF IF YOU WANT TO TEST THIS MODEL
%save nn_search5.mat mean_train sd_train idxavg conversionInfo ...
%    net numfeatureIdx numfeatures miCnst