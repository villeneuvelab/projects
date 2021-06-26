function plotScatterRmse(y,r,t)

    if size(y,2)>size(y,1)
        y = y';
    end
    if size(r,2)>size(r,1)
        r = r';
    end

    figure; hold on;
    mn = floor(0.9*min([y;r])); mx = ceil(1.1*max([y;r]));
    % plot the prediction
    plot(mn:mx,mn:mx,'b--');
    plot(y,r,'r.','MarkerSize',10);
    xlabel('age (years)'); ylabel('predicted age (years)');
    xlim([mn mx]); ylim([mn mx]);
    title(t);

    rmse = sqrt(mean((y' - r').^2));
    text(1.2*mn,0.9*mx,['rmse (unity) = ' num2str(rmse)])
    
    % plot the fit
    p = polyfit(y,r,1);
    py1 = (p(1)*mn) + p(2);
    py2 = (p(1)*mx) + p(2);
    plot([mn mx],[py1 py2],'k','LineWidth',1);
    text(1.2*mn,0.85*mx,['r2 = ' num2str(corr(y,r)^2)]);
    
    legend({'unity';'data';'linear fit'})
end