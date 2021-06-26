
clear
clc

file = '/home/.../.../';  % Load file containing the time series for each participant
mut = -1;
group = 'test';
studygroup = 'all';
ylab = 'age';
[id,y,study,fd,mutation] = load_data(file,mut,group,studygroup,ylab);

if strcmp(group,'test')
    load TrainInfo_extra_new.mat
end

disp('selecting links ...')
if strcmp(group,'train')
    linkInfo = selectConnections('each','each',-1);
else
    linkInfo = TrainInfo.linkInfo;
end

disp('creating feature matrix ...')
FeatureMatInfo = ...
    makeGraphFeatureMat_fromTimeSeries(id,study,y,fd,...
    linkInfo.Pairs,mutation);

%}
if strcmp(group,'train')
    
    disp('doing MR regression on node-wise GraphData ... ')
    [GraphData, MRGraphRegressionInfo] = ...
        regressMRMatrix(FeatureMatInfo.GraphData,[]);
    
    disp('doing community analysis ... ')
    crossModInfo = ...
        getCrossModCorrelations(GraphData,[],linkInfo.Pairs,1,[]);
    
else %test MR regression and communities done based on training data
    
    disp('doing MR regression on node-wise GraphData using training data... ')
    [GraphData, MRGraphRegressionInfo] = ...
        regressMRMatrix(FeatureMatInfo.GraphData, ...
        TrainInfo.MRGraphRegressionInfo);
    
    disp('doing community analysis using training communities... ')
    crossModInfo = ...
        getCrossModCorrelations(GraphData,[],linkInfo.Pairs,1, ...
        TrainInfo.crossModInfo.groupMods);
    
    clear TrainInfo
end


GraphMetrics = [];
for i = 1:size(GraphData,3)
    tic
    disp(['calculating graph metrics for ' int2str(i) ' of ' ...
        int2str(size(GraphData,3))]);
    D = squeeze(GraphData(:,:,i));
    metrics = graph_gamut(D,0.05); %use 5 percent density
    names = fieldnames(metrics);
    gm = zeros(1,length(names));
    for j = 1:length(names)
        eval(['gm(j) = metrics.' names{j} ';']);
    end
    GraphMetrics(i,:) = gm;
    %save tempgraphmetrics.mat GraphMetrics i
    toc
end

% put it all into 1 struct variable
if strcmp(group,'train')
    TrainInfo = [];
    TrainInfo.FeatureMatInfo = FeatureMatInfo;
    TrainInfo.names = names;
    TrainInfo.GraphMetrics = GraphMetrics;
    TrainInfo.linkInfo = linkInfo;
    TrainInfo.MRGraphRegressionInfo = MRGraphRegressionInfo;
    TrainInfo.crossModInfo = crossModInfo;
    %save TrainInfo_extra_new.mat TrainInfo
end

if strcmp(group,'test')
    TestInfo = [];
    TestInfo.FeatureMatInfo = FeatureMatInfo;
    TestInfo.names = names;
    TestInfo.GraphMetrics = GraphMetrics;
    TestInfo.linkInfo = linkInfo;
    TestInfo.MRGraphRegressionInfo = MRGraphRegressionInfo;
    TestInfo.crossModInfo = crossModInfo;
    save TestInfo_extraPAD.mat TestInfo
    %save TestInfo_extra_new.mat TestInfo
end
%}