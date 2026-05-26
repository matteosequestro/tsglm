function f1 = tsglm_plot_individual_pars(modelout)

parnames_to_plot = modelout.par_names;
parnames_to_plot = cellfun(@(x) replace(x, '_', ' '), parnames_to_plot, 'UniformOutput', false);

if nargin < 3
    use_par_ci = 0;
end


n_pars = size(modelout.pars_series, 3);

nrows_tiles = ceil(sqrt(n_pars));
ncols_tiles = ceil(n_pars / nrows_tiles);

% Prepare the layout (keep it tight)
f1 = figure;
tld             = tiledlayout(nrows_tiles, ncols_tiles);
tld.TileSpacing = 'compact';
tld.Padding     = 'compact';

% Common labels
ylabel(tld, 'parameter estimate (a.u.)');
xlabel(tld, 'time (s)');

for parn = 1 : n_pars
    indest = modelout.pars_series(:,:,parn);

    tse = std(indest) ./ sqrt(height(indest));

    time = 1 : width(indest);
    nexttile
    plot(indest', 'Color', [0 0 0 .2])
    hold on
    % patch([time, fliplr(time)], [lowerbound, fliplr(upperbound)], ...
    % [.5,.5,.5], 'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility', 'off')
    % plot(groupest, 'Color', [0 0 0], 'LineWidth', 2)
    yline(0, 'r--', 'LineWidth', 1.5)
    xlim([1, length(time)])

    title(parnames_to_plot{parn})




end


end