function plotScatterRmse_breakout(y,r,t,ptrain)

    if size(y,2)>size(y,1)
        y = y';
    end
    if size(r,2)>size(r,1)
        r = r';
    end

    figure; hold on;
    % plot unity
    mn = floor(0.9*min([y;r])); mx = ceil(1.1*max([y;r]));
    plot(mn:mx,mn:mx,'b');
    
    % plot the test set prediction and linear fit 
    plot(y,r,'r.','MarkerSize',10);
    p = polyfit(y,r,1);
    py1 = (p(1)*mn) + p(2);
    py2 = (p(1)*mx) + p(2);
    plot([mn mx],[py1 py2],'r--','LineWidth',1);
    text(1.2*mn,0.85*mx,['r2 = ' num2str(corr(y,r)^2)]);
    
    % plot the training linear fit
    x = mn:mx; fx = (x*ptrain(1))+ptrain(2);
    plot(x,fx,'k--');
    
    
    xlabel('age (years)'); ylabel('predicted age (years)');
    xlim([mn mx]); ylim([mn mx]);
    title(t);
    rmse = sqrt(mean((y' - r').^2));
    text(1.2*mn,0.9*mx,['rmse (to unity) = ' num2str(rmse)])
    legend({'unity';'test data';'test fit';'train fit'})
end