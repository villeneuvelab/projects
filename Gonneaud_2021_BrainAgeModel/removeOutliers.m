function outliers = removeOutliers(X,obsdim,prc)
    
    %flip so observations are on first axis
    if obsdim ~= 1
        X = X';
    end
    Report = zeros(size(X));
    n = size(X,2); %variables
    
    outliers = [];
    for i = 1:n
        
        x = X(:,i);
        finf = find(isinf(x));
        fnan = find(isnan(x));
        outliers = [outliers; finf; fnan];
        
        if ~isempty(prc)
            s = std(x);
            fsu = find(x>(mean(x) + (prc*std(x))));
            fsl = find(x<(mean(x) - (prc*std(x))));
            outliers = [outliers; fsu; fsl];
            Report(fsu,i) = 1;
            Report(fsl,i) = 1;
        end
        
    end
    
    outliers = unique(outliers);
    %disp([int2str(length(outliers)) ' outliers in this data set']);
    
    if ~isempty(outliers)
        disp('OUTLIERS: ');
        for i = 1:length(outliers)
            x = Report(outliers(i),:);
            y = find(x);
            disp(['observation ' int2str(outliers(i)) ' at columns: '])
            disp(y)
        end
        
        %[~,c] = find(Report);
        %figure;% h = histc(c,1:n); bar(h);
        figure; imagesc(Report); xlabel('column'); ylabel('observation');
        title('outlier locations');
        
    end
    
    
end