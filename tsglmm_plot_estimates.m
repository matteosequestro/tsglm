function f1 = tsglmm_plot_estimates(modelout, cfg, layout_dims)


parnames = modelout.pars.parnames;
tslen = size(modelout.pars.estimates, 2);
npar = length(parnames);
pars_series = modelout.full_ind_estimates;
obs_clusters_sum = modelout.obs_clusters_sum;
% Remove the underscores from the parameter names otherwise it's ugly
parnames_to_plot = cellfun(@(x) replace(x, '_', ' '), parnames, 'UniformOutput', false);

% Time vector
time = 1 : tslen;

% Number of columns and rows for the layout
if nargin < 3
    nrows_tiles = ceil(sqrt(npar));
    ncols_tiles = ceil(npar / nrows_tiles);
else
    nrows_tiles = layout_dims(1);
    ncols_tiles = layout_dims(2);
end
% Prepare the layout (keep it tight)
f1              = figure;
tld             = tiledlayout(nrows_tiles, ncols_tiles);
tld.TileSpacing = 'compact';
tld.Padding     = 'compact';

% Common labels
ylabel(tld, 'parameter estimate (a.u.)');
xlabel(tld, 'time (s)');

% Plot each parameter
for par = 1 : npar
    nexttile

    % Compute stats to plot
    this_parameter      = pars_series(:,:,par);
    grand_average       = mean(this_parameter, 1);
    standard_error      = std(this_parameter, [], 1)  ./ sqrt(height(this_parameter));  % standard error for the confidence interval
    ts                  = tinv([0.025  0.975], height(this_parameter)-1);               % t-score for the confidence interval
    lower_bound         = grand_average + ts(1) .* standard_error;                      % lower bound of confidence interval
    upper_bound         = grand_average + ts(2) .* standard_error;                      % upper bound of confidence interval

    % The shade (intervals)
    patch([time, fliplr(time)], [lower_bound, fliplr(upper_bound)], ...
        [0.5 0.5 0.5],'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
    hold on

    % Grand average
    plot(grand_average, 'Color', [0 0 0], 'LineWidth', 2)
    if isfield(cfg, "stim_onset_time")
        xline(cfg.stim_onset_time, '--')
    end

    % Plot the group-level (hyper) parameters if required. Note that
    % they should be the same as the mean of individual estimates, so
    % plotting them is mostly for diagnostic. If they are the 2 time
    % series do not concide something is wrong
    if isfield(cfg, "plot_group_pars")
        if cfg.plot_group_pars
            plot(modelout.pars.estimates(par,:), 'r--', 'LineWidth', 1.5);
        end
    end

    % Put a title with the parameter name
    title(parnames_to_plot{par})

    % Draw an horizontal line on y = 0 as a reference
    yline(0, '--', 'HandleVisibility', 'off')

    % Draw a rectangle for each cluster, the color depends on the
    % significance
    this_par_clusters = modelout.obs_clusters_sum(strcmp(obs_clusters_sum.pred, parnames{par}),:);
    for this_cluster = 1 : height(this_par_clusters)
        hold on
        clFirst = this_par_clusters.first(this_cluster);
        clLength = this_par_clusters.length(this_cluster)-1;
        ylims = get(gca, 'Ylim');
        if this_par_clusters.prob(this_cluster) >= modelout.alpha
            if cfg.plot_non_significant_clusters
                rectColor = [.5 .5 .5];
                rectangle('Position', [clFirst, ylims(1), clLength, sum(abs(ylims))], ...
                    'FaceColor', [rectColor .10], 'EdgeColor',[rectColor .1], 'FaceAlpha', .2)
            end
        else
            rectColor =[0.3569 0.6078 0.8353];
            rectangle('Position', [clFirst, ylims(1), clLength, sum(abs(ylims))], ...
                'FaceColor', [rectColor .10], 'EdgeColor',[rectColor .1], 'FaceAlpha', .2)
        end


    end

    % Cut the x axis to exclude parts of signal you don't want
    xlim([1, max(time)]);

    % Adjust the time axis
    ax = gca; % Get current axis handle
    if isfield(cfg, "stim_onset_time")
        ax.XTickLabel = (ax.XTick - cfg.stim_onset_time) / cfg.fs; % Divide by sample frequency
    else
        ax.XTickLabel = ax.XTick / cfg.fs; % Divide by sample frequency
    end

    % Draw a box cause I like it
    box on

end % for par = 1 : par


end