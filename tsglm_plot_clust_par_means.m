function [out, f1] = tsglmm_plot_clust_par_means(data, version, verbose)
% plots group-level parameter means with SEs, Bayes factors (BF10),
% and combined dual-axis view.

if nargin < 3
    verbose = 1;
end

parnames = data.Properties.VariableNames;    % get parameter names
parnames = cellfun(@(x) replace(x, "_", " "), parnames, 'UniformOutput', false);

data = table2array(data);                    % convert table to numeric array

means = mean(data);                          % compute mean across participants
ses   = std(data) ./ sqrt(height(data));     % compute standard error of the mean

%%% Compute Bayes factors for each parameter (vs zero)
BF10 = NaN(width(data), 1);
for ii = 1:width(data)
    BF10(ii) = bf.ttest(data(:,ii), 0);      % one-sample Bayesian t-test
    if verbose
        fprintf('%s: BF10 = %.2f, BF01 = %.2f \n', ...
            parnames{ii}, BF10(ii), 1/BF10(ii));
    end
end

%%% Plot Bayes factors and parameter means

if strcmp(version, "long")
    f1 = figure;

    tiledlayout(1,3)                             % 3 panels side by side

    % -------------------------------------------------------------------------
    % (1) Individual averages for each cluster
    % -------------------------------------------------------------------------
    nexttile
    hold on

    b = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', ...
        'FaceAlpha', 0.3);
    errorbar(1:numel(means), means, ses, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);
    for i = 1:width(data)
        x = i + 0.1 + 0.1*(rand(size(data,1),1)-0.5);  % add small horizontal jitter
        scatter(x, data(:,i), 15, [0.5 0.5 0.5], 'filled', ...
            'MarkerFaceAlpha', 0.3);             % individual dots
    end

    xlim([0.5 width(data)+.5])
    set(gca, 'XTick', 1:width(data), 'XTickLabels', parnames)
    box on
    yline(0, '--')
    ylabel('parameter estimate (a.u.)')
    title('By-subject cluster averages')

    % -------------------------------------------------------------------------
    % (2) Bayes factors (log10 scale)
    % -------------------------------------------------------------------------
    nexttile
    hold on
    scatter(1:width(data), log10(BF10), 50, 'r', 'filled')
    yline(log10(3),  '--k', 'BF10 = 3 (moderate)', 'LabelVerticalAlignment','bottom')
    yline(log10(10), '--k', 'BF10 = 10 (strong)',  'LabelVerticalAlignment','bottom')

    set(gca, 'XTick', 1:width(data), 'XTickLabels', parnames)
    ylabel('log_{10}(BF_{10})')
    xlim([0.5 width(data)+0.5])
    box on
    title('Evidence for parameter ≠ 0')

    % -------------------------------------------------------------------------
    % (3) Combined plot: means ± SE and log10(BF10) with dual y-axis
    % -------------------------------------------------------------------------
    nexttile
    hold on
elseif (strcmp(version, "noplot"))
    %%% Output structure
    out.parnames = parnames;
    out.BF10     = BF10;
    out.BF01     = 1 ./ BF10;
    out.means    = means;
    out.ses      = ses;

    f1 = nan;
else
    f1 = figure;

    hold on

    for i = 1:width(data)
        x = i + 0.1*(rand(size(data,1),1)-0.5);  % add small horizontal jitter
        scatter(x, data(:,i), 15, [0.5 0.5 0.5], 'filled', ...
            'MarkerFaceAlpha', 0.3);             % individual dots
    end
    errorbar(1:width(data), means, ses, 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', [1 0 0], 'Color', [1 0 0], ...
        'LineWidth', 1.2, 'CapSize', 0);
    ylabel('mean estimate in cluster (a.u.)')
    yline(0, '--')

    set(gca, 'XTick', 1:width(data), 'XTickLabels', parnames)
    xlim([0.5 width(data)+.5])
    box on

    %%% Output structure
    out.parnames = parnames;
    out.BF10     = BF10;
    out.BF01     = 1 ./ BF10;
    out.means    = means;
    out.ses      = ses;
end
end
