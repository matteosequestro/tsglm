function f1 = nice_simple_corrplot(x,y, xlab, ylab, tcolor) 
    
% Set defaults
if ~exist("xlab", 'var'); xlab = ""; end
if ~exist("ylab", 'var'); ylab = ""; end
if ~exist("tcolor", 'var'); tcolor = [0 .5 .5]; end

% Fit linear model
mdl = fitlm(x, y);

xq = linspace(min(x), max(x), 100)'; % sorted x values for smooth line

[yhat, yci] = predict(mdl, xq);


% figure;
if nargout > 0
    f1 = figure;
end

plot(xq, yhat, 'Color', tcolor, 'LineWidth', 2)
hold on

fill([xq; flipud(xq)], [yci(:,1); flipud(yci(:,2))], tcolor, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xq, yhat, 'Color',tcolor, 'LineWidth', 2)

% Scatter plot
scatter(x, y, 'filled', 'MarkerFaceColor', tcolor, ...
        'MarkerFaceAlpha', 0.5, 'SizeData', 30); 
xlabel(xlab)
ylabel(ylab)

box on;


end

