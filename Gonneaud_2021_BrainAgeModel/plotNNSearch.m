function plotNNSearch(x,numfeatures,NETS,t)

    figure; imagesc(x'); 
    ystr = {};
    for i = 1:size(NETS,2)
        ystr{i} = num2str(NETS{1,i});
    end
    set(gca,'XTick',1:size(x,1)); set(gca,'XTickLabels',numfeatures);
    set(gca,'YTick',1:size(x,2)); set(gca,'YTickLabels',ystr);
    title(t); ylabel('network'); xlabel('# features')
    
end