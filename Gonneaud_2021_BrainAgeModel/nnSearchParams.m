function [PERFORMANCE, NETS, numfeatures] = ...
    nnSearchParams(DataRank,DataTestRank,y,yTest)

    % Now take the ranked features and test them in multiple neural network
    % architectures
    numfeatures = 5:5:size(DataRank,2); % number of inputs (features)
    u1 = 5:5:15;  % number of units in the first hidden layer
    u2 = [0 2 5]; % number of units in the second hidden layer

    PERFORMANCE = [];%zeros(length(u1)*length(u2)-1);
    NETS = {};

    randiter = 3; % each network is tested 3 times, and performance is averaged
    for i = 1:length(numfeatures)
        Data0 = DataRank(:,1:numfeatures(i));
        DataTest0 = DataTestRank(:,1:numfeatures(i));

        ct = 0;
        for j = 1:length(u1)

            for jj = 1:length(u2)
                units = [u1(j) u2(jj)];

                if units(2) < units(1)
                    if units(2)==0 
                        units = units(1);
                    end

                    tic

                    rp = randperm(size(DataRank,2));
                    DataNull0 = DataRank(:,rp(1:numfeatures(i)));

                    disp(['training on ' int2str(numfeatures(i)) ...
                        ' features with network ' int2str(units)]);

                    ct = ct+1;
                    NETS{i,ct} = units;

                    rmse_train = []; r2_train = [];
                    rmse_test = []; r2_test = [];
                    rmse_null = []; r2_null = [];
                    for jjj = 1:randiter

                        [net,ynntrain] = tryNN(Data0,y,units,0);
                        rmse_train(jjj) = sqrt(mean((ynntrain' - y).^2));
                        r2_train(jjj) = corr(ynntrain',y)^2;

                        ynntest = net(DataTest0');
                        rmse_test(jjj) = sqrt(mean((ynntest' - yTest).^2));
                        r2_test(jjj) = corr(ynntest',yTest)^2;

                        [~,ynnnull] = tryNN(DataNull0,y,units,0);
                        rmse_null(jjj) = sqrt(mean((ynnnull' - y).^2));
                        r2_null(jjj) = corr(ynnnull',y)^2;

                    end

                    PERFORMANCE.rmse_train(i,ct) = mean(rmse_train);
                    PERFORMANCE.r2_train(i,ct) = mean(r2_train);
                    PERFORMANCE.rmse_test(i,ct) = mean(rmse_test);
                    PERFORMANCE.r2_test(i,ct) = mean(r2_test);
                    PERFORMANCE.rmse_null(i,ct) = mean(rmse_null);
                    PERFORMANCE.r2_null(i,ct) = mean(r2_null);

                    toc
                end
            end
        end
    end
end