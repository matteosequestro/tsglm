function f1 = plot_within_error(fin, errortype, colors, legendlabs)



% Provide default color if not provided as argument. Give an error if
% they are not enough to cover all conditions
if ~exist('colors', 'var')
    colors = ['#FF0000'; '#00FF00'; '#0000FF'; '#FF00FF'; '#00FFFF'; ...
        '#800000'; '#008000'; '#000080'; '#FFA500'; '#A52A2A'];
    if height(colors) < height(fin)
        error(['You need to provide ', num2str(height(fin)), ' colors, but you only have ', num2str(height(colors)) ' by default. You need to provide colors manually' ]);
    end
else
    if height(colors) < height(fin)
        error(['You need to provide ', num2str(height(fin)), ' colors, but you only have ', num2str(height(colors))  ]);
    end
end

% I couldn't make the function I use for plotting CIs to accept
% hexadecimal colors, so this is converting them in RGB
colors_rgb = cell(height(colors), 1);
for col = 1 : height(colors)
    colors_rgb(col)  = {sscanf(colors(col, 2:end),'%2x%2x%2x',[1 3])/255};
end

% Big Line (mean)
mean_line = fin.mean_series;

% Define the patches depending on which error you want
if ~exist('errortype', 'var'); errortype = "se"; end % set standard error as default error type
if ~exist('wantpatch', 'var'); wantpatch = 1; end % set standard error as default error type

if strcmp(errortype, "se")
    up_lim = fin.se_up;
    low_lim = fin.se_low;
elseif strcmp(errortype, "ci")
    up_lim = fin.ci_up;
    low_lim = fin.ci_low;
end

% Plot each condition
if nargout ~= 0
    f1 = figure;
end
time = 1:width(fin.mean_series);
for rr = 1  : height(fin)
    % Define legend name for this condition
    dispname= '';
    for nn = 1 : (width(fin) - 5)
        dispname = [dispname, [ fin.Properties.VariableNames{nn} ': ' num2str(fin{rr,nn}) ,' | '] ];
    end
    dispname((end-2):end) = [];

    % Plot this condition
    if wantpatch
        patch([time, fliplr(time)], [low_lim(rr,:), fliplr(up_lim(rr,:))], ...
            colors_rgb{rr},'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
    end


    hold on
    plot(mean_line(rr,:), 'Color', colors_rgb{rr}, 'LineWidth', 1, 'DisplayName',dispname)
end
xlim([1, max(time)]);
if nargin < 4
    legend('Location', 'north')
else
    legend(legendlabs, 'Location', 'north')
end

box on


end