%%% compute averages and within errors (se and ci) over time-series, then
%%% plot them

function [fin, exp_ind_means, f1] = within_error(set, series, pID, conditions, ci_lims, want_plot, wantpatch, errortype, colors,legendlabs)
%%% ========================================================================================================================
%%% set         : full dataset having at least a column for ID, x columns for the condition(s) you want to average over, and one column with the time series in a cell. Each row is a trial
%%% series      : the name of column with the time series
%%% pID         : the nale of column with the IDs (or the clustering variable in general)
%%% conditions  : string array with the name of the names of the conditions you want to average over (e.g., ["condition", "valence"])
%%% want_plot   : 0 if you don't want a plot of the series, 1 if you do

%%% errortype   : "se" for standard error, "ci" for confidence interval. This
%%% is just for plotting, in the output table you will have both

%%% colors      : vertical array with hexadecimal colors (one for each condition
%%% you want to average over). This is just for plotting, it should look
%%% like this: colors = ['#C0C0C0';'#808080'; '#000000';'#FFA500'; '#A52A2A';'#800000'];

%%% ci_lims     : limits for the confidence interval. default is ci_lims = [0.025, 0.975] for a 95% CI


%%% EXAMPLE OF USE:
%%% within_error(out, "cutInterpHrLinearBC", "id", ["condition", "valence"], [], 1, 1);
%%% or
%%% within_error(out, "cutInterpHrLinearBC", "id", "valence", [], 1, 1);
%%% ========================================================================================================================


% Set Defaults ------------------------------------------------------------
if ~exist("ci_lims", 'var') || isempty(ci_lims); ci_lims = [0.025, 0.975]; end
if ~exist("want_plot", 'var') || isempty(want_plot); want_plot = 1; end
if ~exist("errortype", 'var') || isempty(errortype); errortype = "ci"; end
if ~exist("legendlabs", 'var') || isempty(errortype); errortype = "ci"; end

idvar= pID;


% find the columns for conditions -----------------------------------------
series = set{:, series};
pID = set{:,pID};
conditions_tot = zeros(height(set), length(conditions));
for ll = 1 : length(conditions)
    if iscell(set{:, conditions(ll)})
        conditions_tot(:, ll) = cell2mat(set{:, conditions(ll)});
    else
        conditions_tot(:, ll) = set{:, conditions(ll)};
    end
end

% averages by-participant -------------------------------------------------
upID = unique(pID);
ind_means = zeros(length(upID), length(series{1}));
for pp = 1 : length(upID)
    this_pp = series(strcmp(pID, upID(pp)), :);
    mean_pp = nanmean(cell2mat(this_pp), 1);
    ind_means(pp, :) = mean_pp;
end

% grand-averag%e ------------------------------------------------------------
grand_mean = nanmean(ind_means, 1);

% adjustment factor (grand average - individual average) ------------------
adj_factor  = grand_mean - ind_means;
tadj_factor = series;
for ii = 1 : height(series)
    which_pid =strcmp(upID, pID(ii));
    which_adj = adj_factor(which_pid, :);
    tadj_factor{ii} = which_adj ;
end

% adjusted values for each variable (same the mean, adjusted error) -------
adjusted_series = cell2mat(series) + cell2mat(tadj_factor);

%% Means for each condition --------------------------------------------------------------------------------------------------

% Extract possible combinations of conditions
conds_cat = struct([]);
for kk = 1 : width(conditions_tot)
    tcond = [conditions_tot(:, kk)];
    conds_cat(kk).values = [unique(tcond)];
end

conds_combinations = cell(size(conds_cat));
[conds_combinations{:}] = ndgrid(conds_cat.values);
conds_combinations = array2table(cell2mat(cellfun(@(x) x(:), conds_combinations, 'UniformOutput', false)));
conds_combinations.Properties.VariableNames = conditions;


% Mean for each combination
exp_ind_means = [];
for pp = 1 : length(upID)
    this_pp = adjusted_series(strcmp(set.(idvar), upID(pp)), :);
    conditions_pp = conditions_tot(strcmp(set.(idvar), upID(pp)), :);
    for rr = 1 : height(conds_combinations)
        export = cell(1, width(conds_combinations) + 2);

        try export(1, 1) = upID(pp);
        catch; export{1, 1} = upID(pp); end

        export(1, 2 : (width(conds_combinations)+1)) = table2cell(conds_combinations(rr, :));

        this_comb  = {this_pp(ismember(conditions_pp, conds_combinations{rr, :}, 'rows'), :)};
        export(1, end) = {nanmean(cell2mat(this_comb), 1)};

        exp_ind_means = [exp_ind_means; export];
    end
end
exp_ind_means = cell2table(exp_ind_means);
exp_ind_means.Properties.VariableNames = [{'id'} conds_combinations.Properties.VariableNames {'mean_series'}];



%% Compute within errors
fin = cell(height(conds_combinations), (5 + width(conds_combinations)));
nconds =  height(conds_combinations);
for rr = 1 : nconds

    this_cond        = exp_ind_means(ismember(exp_ind_means(:, conds_combinations.Properties.VariableNames), conds_combinations(rr, :)), :);
    this_cond_mean   = nanmean(this_cond.mean_series, 1);


    this_cond_se     = nanstd(this_cond.mean_series) / sqrt(height(this_cond));
    this_cond_se     = this_cond_se * sqrt(nconds / (nconds - 1));

    this_cond_se_low = this_cond_mean - this_cond_se;
    this_cond_se_up  = this_cond_mean + this_cond_se;

    ts               = tinv([ci_lims(1)  ci_lims(2)], height(this_cond)-1);      % T-Score
    this_cond_ci_low = this_cond_mean + ts(1) * this_cond_se;
    this_cond_ci_up  = this_cond_mean + ts(2) * this_cond_se;

    fin(rr, :)       = table2cell([this_cond(1, 2:end-1) {this_cond_mean}, {this_cond_se_up}, {this_cond_se_low}, {this_cond_ci_up}, {this_cond_ci_low}]);
end
fin = cell2table(fin);
fin.Properties.VariableNames = [conds_combinations.Properties.VariableNames, {'mean_series'}, {'se_up'}, {'se_low'},{'ci_up'},{'ci_low'}];



%% Plot (if you want)

if want_plot
    if nargout == 0
        plot_within_error(fin, errortype, colors, legendlabs);
    else
        try
            f1  =  plot_within_error(fin, errortype, colors, legendlabs);
        catch
            f1  = plot_within_error(fin, errortype);
        end
    end
else
    f1 =  [];
end


end