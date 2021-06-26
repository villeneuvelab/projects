% returns index of features which seem to be more important, ignoring a
% handful of the others which tend to screw things up

function keepmetric = filterFeatures(names)

    keepmetric = []; ct = 0;
    for i = 1:length(names)
        if isempty([strfind(names{i},'iqr'), ...
                strfind(names{i},'25'), strfind(names{i},'75')])
            ct = ct+1;
            keepmetric(ct) = i;
        end
    end
    
end