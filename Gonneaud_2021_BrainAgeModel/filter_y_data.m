function [y,finidx,idfin] = filter_y_data(file,id,study,ylab)

    [num_data,txt_data,raw_data] = xlsread(file);
    vlabels = txt_data(1,:);
    
    %id_xls = convertId(num_data(:,1),txt_data(2:end,1));
    id_xls = raw_data(2:end,1);
    %study_xls = txt_data(2:end,12);
    study_xls = raw_data(2:end,12);
    %y_xls = num_data(:,find(strcmp(vlabels,ylab)));
    y_xls = cell2mat(raw_data(2:end,find(strcmp(vlabels,ylab))));
    
    m = length(id_xls);
    finidx = [];
    y = [];
    idfin = {};
    ct = 0;
    for i = 1:m
        idx = find(strcmp(string(id),string(id_xls{i})) .* strcmp(study,study_xls{i}));
        if ~isempty(idx)
            ct = ct+1;
            finidx(ct) = idx(1);  %THIS SHOULD NOT HAVE MORE THAN 1 INDEX
            y(ct) = y_xls(i);
            idfin{ct} = id_xls{i};
        end
    end
    
    
    finidx = finidx(~isnan(y))';
    idfin = idfin(~isnan(y));
    y = y(~isnan(y))';
end