function [dataN,m,s] = standardizeData(Data,m,s)

    l = size(Data,1);
    if isempty(m)
        m = mean(Data);
    end
    if isempty(s)
        s = std(Data);
    end
    dataN = (Data-repmat(m,l,1))./repmat(s,l,1);
    
end