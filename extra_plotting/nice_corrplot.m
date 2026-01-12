function f1 = nice_corrplot(x,y, xlab, ylab, titles)

if ~exist("xlab", "var"); xlab = "x";end
if ~exist("ylab", "var"); ylab = "y";end
if ~exist("titles", "var"); titles= "";end

if width(x) == 1 && width(y) == 1
    ccc     = round(concordance_correlation_coefficient(x,y), 2);           % Concordance Correlation Coefficient
    cc      = round(corr(x,y), 2);                                          % Simple pearson correlation
    boot_ci = corr_boot(x,y);                                               % Bootstrap confident interval for the correlation
    
    % Plot the rest
    mdl = fitlm(x,y);

    close all
    f1 = figure;
    scatter(x,y, 'o', 'filled', 'MarkerFaceColor', [255 178 102]/255, ...
        'MarkerFaceAlpha', 0.8, 'SizeData', 30)
    hold on

    xq = linspace(min(x), max(x), 100)'; % sorted x values for smooth line

    [yhat, yci] = predict(mdl, xq);

    fill([xq; flipud(xq)], [yci(:,1); flipud(yci(:,2))], [255 178 102]/255, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(xq, yhat, 'Color', [255 178 102]/255, 'LineWidth', 2)
    % h = plot(mdl, 'Color', [255 178 102]/255);
    xy = [min([xlim, ylim]) , max([xlim, ylim])];
    xlim(xy)
    ylim(xy)
    plot(xy(1):.01:xy(2), xy(1):.01:xy(2), '--', 'Color', [.5 .5 .5])
    b = gca; legend(b,'off');
    xlabel(xlab)
    ylabel(ylab)
    xL=xy;
    yL=xy;

    % centerY = (abs(yL(2)) - abs(yL(1))) / 2;
    poles = [yL(1),yL(2)];
    distance = diff(poles);
    distances = (min(poles) :(abs(distance) * .1) : max(poles));

    Y1 = quantile(distances, .2) ;
    Y2 = quantile(distances, .1) ;
    text(0.95*xL(2),Y1, sprintf('r=%.2f [%.2f, %.2f]', cc, boot_ci(1), boot_ci(2)),'HorizontalAlignment','right','VerticalAlignment','bottom')
    text(0.95*xL(2), Y2,['ccc=', num2str(ccc)],'HorizontalAlignment','right','VerticalAlignment','bottom')
    box on
    title('')
else
    ncol = 3;
    nrow = ceil(width(x) / ncol);
    f1 = tiledlayout(nrow, ncol, 'TileSpacing','compact');
    for pp = 1 : width(x)
        nexttile
        ccc = round(concordance_correlation_coefficient(x(:,pp), y(:,pp)), 2);
        cc = round(corr(x(:,pp), y(:,pp)), 2);
        boot_ci = corr_boot(x(:,pp),y(:,pp));                                               % Bootstrap confident interval for the correlation


        mdl = fitlm(x(:,pp), y(:,pp));

        scatter(x(:,pp), y(:,pp), 'o', 'filled', 'MarkerFaceColor', [255 178 102]/255, 'MarkerFaceAlpha', 0.8, 'SizeData', 15)
        hold on

        xq = linspace(min(x(:, pp)), max(x(:, pp)), 100)'; % sorted x values for smooth line

        [yhat, yci] = predict(mdl, xq);

        fill([xq; flipud(xq)], [yci(:,1); flipud(yci(:,2))], [255 178 102]/255, ...
            'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(xq, yhat, 'Color', [255 178 102]/255, 'LineWidth', 2)
        % h = plot(mdl, 'Color', [255 178 102]/255);
        xy = [min([xlim, ylim]) , max([xlim, ylim])];
        xlim(xy)
        ylim(xy)
        plot(xy(1):.01:xy(2), xy(1):.01:xy(2), '--', 'Color', [.5 .5 .5])
        b = gca; legend(b,'off');
        xlabel(xlab)
        ylabel(ylab)
        xL=xy;
        yL=xy;

        % centerY = (abs(yL(2)) - abs(yL(1))) / 2;
        poles = [yL(1),yL(2)];
        distance = diff(poles);
        distances = (min(poles) :(abs(distance) * .1) : max(poles));

        Y1 = quantile(distances, .2) ;
        Y2 = quantile(distances, .1) ;
        text(0.95*xL(2), Y1,sprintf('r=%.2f [%.2f, %.2f]\n', cc, boot_ci(1), boot_ci(2)),'HorizontalAlignment','right','VerticalAlignment','bottom')
        % text(0.95*xL(2), Y1,sprintf('r=%.2f', cc),'HorizontalAlignment','right','VerticalAlignment','bottom')
        text(0.95*xL(2), Y2,['ccc=', num2str(ccc)],'HorizontalAlignment','right','VerticalAlignment','bottom')
        box on

        if exist("titles", "var")
            title(titles(pp))
        else
            title('')
        end

    end


end

end