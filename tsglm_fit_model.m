%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to compute parameters with GLME for a time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelout = tsglm_fit_model(data, formula, cfg)
% --------------------------------------------------------------------------
% Computes linear mixed-effects models (GLME or LME) across a time series,
% optionally in parallel, and extracts parameter estimates, R-squared,
% information criteria, residuals, and participant-level estimates.
%
% INPUTS:
%   data                  : table containing at least:
%                             - column for ID
%                             - columns for each predictor in the model
%                             - column with time series data in a cell
%                               array. Each cell contains the time series in
%                               an horizontal vector (e.g., if you have
%                               800 samples for each trial, each row should
%                               look like {[1x800 double]}
%                           Each row of the data table corresponds to a trial.
%   formula               : string, model formula in Wilkinson notation,
%                           e.g., 'y ~ x1 * x2 + (1 + x1 * x2 | id)'. At
%                           the moment (September 17th 2025) you need full
%                           random slopes structure for it to work.
%   cfg                   : struct with optional fields:
%       .want_parallel_fit: logical (1/0) to fit models in parallel
%       .nworkers         : number of workers for parallel computation
%       .verbose_fit      : logical (1/0) print progress messages
%       .want_diagnostic  : logical (1/0) compute predictions and residuals
%       .glm_likelihood   : likelihood type for fitglme (default 'Gaussian')
%
% OUTPUTS:
%   modelout              : struct with fields:
%       .pars              : parameter estimates, CI, t-statistics
%       .criteria          : log-likelihood, AIC, BIC, deviance
%       .rsqrd             : ordinary and adjusted R-squared
%       .full_ind_estimates: participant-level parameter estimates
%       .predictions       : predicted and observed averages
%       .residuals         : residuals and standardized residuals
%
% --------------------------------------------------------------------------

%% Set defaults
% Unpack method variables
verbose_fit       = get_or_default(cfg, 'verbose_fit', 0);
want_diagnostic   = get_or_default(cfg, 'want_diagnostic', 1);
want_parallel_fit = get_or_default(cfg, 'want_parallel_fit', 0);
glm_likelihood    = get_or_default(cfg, 'glm_likelihood', 'Gaussian');
nworkers          = get_or_default(cfg, 'nworkers', []);  % empty = default pool



% Parallel pool setup
if want_parallel_fit
    if ~isempty(nworkers)
        maxWorkers = parcluster('local').NumWorkers;
        if nworkers > maxWorkers
            error("Requested %d workers, but only %d are available in the local cluster profile.", ...
                nworkers, maxWorkers);
        end
    end

    pool = gcp('nocreate');  % get existing pool if any
    if isempty(pool)
        if isempty(nworkers)
            parpool;  % use default max workers
        else
            parpool(nworkers);
        end
    elseif ~isempty(nworkers) && pool.NumWorkers ~= nworkers
        delete(pool);
        parpool(nworkers);
    end
end

%%% Define likelihood function
% You could just always use figlme with 'Distribution', 'Gaussian' when
% needed, but fitlme is faster and preferrable
% Fit the model and check timing
if strcmp(glm_likelihood, 'Gaussian')
    fitfun = @(x,y) fitlm(x, y);
else
    fitfun = @(x,y) fitglme(x, y, 'Distribution');
end

%%% Extract dependent variable
yvar = strtrim(formula(1:find(formula=='~')-1));

%%% Find the time series variable (which is not necessarily the dependent
%%% variable if, e.g., you are predicting choices by signal)
fixed_eff_formula   = formula(find(formula=='~')+1 : end);
fixed_eff_formula   = erase(fixed_eff_formula, ["+", "*"]);
predictors          = strsplit(fixed_eff_formula, ' ');
predictors(cellfun(@(x) isempty(x) || strcmp(x, '1') || contains(x,':'), predictors)) = [];

all_vars = [predictors, {yvar}];


all_tsvars = all_vars(cellfun(@(x) isa(data.(x), 'cell') && size(cell2mat(data.(x)), 2) > 1, all_vars ));
if isempty(all_tsvars)
    error(['Model formula needs to contain at least one time series variable. Current formula is: ', formula])
elseif isscalar(all_tsvars)
    tsvar = all_tsvars{1};
elseif length(all_tsvars) > 1
    error(['your formula contains ' num2str(length(all_tsvars)), ' time series variable, but currently tsglmm_fit_model supports only one. If you really need more you can ask me to hurry up in implementing this feature at matteo.sequestro@gmail.com'])
end

%%% Time series length
tslen = length(data.(tsvar){1});

%%% Backup of the dependent variable (needed later)
data.ts2 = data.(tsvar);

%% Fit one time point to estimate runtime and preallocate arrays
% Extract a single sample for the dependent variable
tmpset      = data;
tmpset.(tsvar) = cellfun(@(x) x(1), tmpset.(tsvar));

% fit one model
tic
tmp_rm = fitfun(tmpset, formula);
fake_length = toc;

% find id variable
idvar = tmp_rm.PredictorNames{~ismember(tmp_rm.PredictorNames, all_vars)};

% build minimal dataset
min_data = data(:, [all_vars, idvar, "ts2"]);

% Display runtime for a single model and estimated overall runtime
expected_length = fake_length * tslen;
if want_parallel_fit
    pool = gcp('nocreate');
    if ~isempty(pool)
        expected_eff = 0.6; % estimated parallel efficiency (by "estimated" I mean that I ran the function a bunch of times and this worked well enough)
        expected_length = expected_length / (pool.NumWorkers * expected_eff);
    end
end

fprintf('One model took %.2f seconds. Rough estimate of total time: %.2f seconds (~%.2f minutes). Adjust expectations accordingly.\n', ...
    fake_length, expected_length, expected_length/60);

% create a tmp folder
rndname = sprintf('tmp_%06d', randi(1e6));
tmpout = fullfile(cd, rndname);
mkdir(tmpout);


%% ============================== fit all other time points ==============================
% the code for parallel and serial loop is the same - I should probably do
% a function and repeat it, but the first time I tried for a few minutes
% MATLAB didn't like and was giving weird errors so I gave up. Cheers.
if want_parallel_fit
    parfor tt = 1:tslen
        if verbose_fit; fprintf("parallel modeling point: %d/%d\n", tt, tslen); end
        % get the time point to fit
        set = min_data;
        set.(tsvar) = cellfun(@(x) x(tt), min_data.ts2);
        set.ts2 = [];

        % fit the model
        rm = fitfun(set, formula);

        % make a name for the tmp file
        printmod = ['%0', num2str(length(num2str(tslen))), 'd'];
        outnum   = sprintf(printmod, tt);
        outname  = tmpout + "\tmpfit" + outnum + ".mat";

        % join what to export and export
        output     = struct();
        output.rm  = rm;

        save(outname,'-fromstruct',  output)
    end
else
    for tt = 1:tslen
        if verbose_fit; fprintf("parallel modeling point: %d/%d\n", tt, tslen); end
        % get the time point to fit
        set = min_data;
        set.(tsvar) = cellfun(@(x) x(tt), min_data.ts2);
        set.ts2 = [];

        % fit the model
        rm = fitfun(set, formula);

        % make a name for the tmp file
        printmod = ['%0', num2str(length(num2str(tslen))), 'd'];
        outnum   = sprintf(printmod, tt);
        outname  = tmpout + "\tmpfit" + outnum + ".mat";

        % join what to export and export
        output     = struct();
        output.rm  = rm;

        save(outname,'-fromstruct',  output)

    end
end



%% join time points
tmpfiles = dir(fullfile(tmpout, 'tmpfit*.mat'));

npars    = length(tmp_rm.CoefficientNames);          % number of parameters (or coefficients if you prefer)
parnames = tmp_rm.CoefficientNames;                  % parameters' names
parnames = sort(parnames);

idvar    = tmp_rm.Formula.GroupingVariableNames{1};  % the name of the column containing participants' ids
nids     = length(unique(data{:,idvar}));            % number of participants


% Compute VIF values (unless it's an only intercept model
if npars > 1
    modelout.vifs   = tsglmm_compute_vif(tmp_rm);                  % Variance Inflaction Factor (VIF)
else
    modelout.vifs = NaN;
end

%%% Preallocate arrays
loglik      = zeros(1, tslen);
AIC         = zeros(1, tslen);
BIC         = zeros(1, tslen);
deviance    = zeros(1, tslen);
ordinary    = zeros(1, tslen);
adjusted    = zeros(1, tslen);


%%% Preallocate other arrays
estimates   = zeros(npars, tslen);              % parameter estimate for each group level parameter and time point
lower_ci    = estimates;                        % lower boundary of the 95% CI for each group level parameter and sample point
upper_ci    = estimates;                        % upper boundary of the 95% CI for each group level parameter and sample point
t_stat      = estimates;                        % t statistics for each group level parameter and time point
full_ind_estimates = NaN(nids, tslen, npars);   % N participants X t samples X V parameters matrix
predictions = NaN(nids, tslen);                 % estimated value of the dependent variable for each time point
averages    = NaN(nids, tslen);                 % average of the dependent variable for each participant
residuals   = NaN(height(data), tslen);         % model residuals for each sample point

parfor tt = 1 : tslen
    if verbose_fit; fprintf("processing time point: %d/%d\n", tt, tslen); end

    tsample = load(fullfile(tmpfiles(tt).folder, tmpfiles(tt).name));
    rm = tsample.rm;
    set = rm.Variables;

    % Participant averages
    averages(:,tt) = groupsummary(set.(yvar), table2array(set(:,idvar)), 'mean');
    
    % Parameter names in this particular fit
    parnames2 = rm.CoefficientNames;                  % parameters' names

    % Sometimes the predictor names have different arrangements in
    % different fits, so fix it (terrible and inelegant solution but hey it
    % works)
    missing = ~ismember(parnames2, parnames);
    to_change = parnames2(missing);
    to_change_later = [];
    for ii = 1 : length(to_change)
        % isolate the variable to process
        this_var = to_change{ii};


        colons = find(this_var == ':');
        colons = [1, colons, length(this_var)];
        vars = [];
        for tc = 2 : length(colons)
            tv = this_var(colons(tc-1):colons(tc));
            tv = replace(tv, ':', '');
            vars = [vars; change_text(tv)];
        end
        
        % find the original one
        which_one = ones(1, length(parnames));
        for pr = 1 : height(vars)
            which_one = which_one .* contains(parnames, vars(pr));
        end

        % exclude the possibility that a 3way predictor is kept despite the
        % variable being 2 ways (another ugly solution but hey it works)
        if sum(which_one) ~=1
            hom_many_colons = cellfun(@(x) sum(x == ':'), parnames);
            true_colons_n = length(colons) -2;
            which_one(hom_many_colons ~= true_colons_n) =0;
        end

        % save the ones to change cause you'll do the same later
        to_change_later = [to_change_later; [parnames2(strcmp(parnames2, to_change(ii))),parnames(logical(which_one)) ]];

        % finally change it
        parnames2(strcmp(parnames2, to_change(ii))) = parnames(logical(which_one));
    end
    
    % sort the names so it matches the original one (already sorted
    parnames2 = sort(parnames2);

    % in fact, if eventually this set is not the same as the original one,
    % throw an error
    if ~all(strcmp(parnames, parnames2))
        error("different fits have different parameter names")
    end

    % export diagnostics if necessary
    if want_diagnostic
        yhat = predict(rm, set);
        predictions(:,tt) = groupsummary(yhat, table2array(set(:,idvar)), 'mean');
        residuals(:, tt) = set.(yvar) - yhat;
    end

    % Extract fixed and random effects
    [~, ~, festats] = fixedEffects(rm);
    [~, ~, blupstats] = randomEffects(rm);
    
    festats = sortrows(festats, "Name");

    % change the parameter names to match the original one
    for pc = 1 : height(to_change_later)
        this_change = to_change_later(pc,:);
        festats.Name(strcmp(festats.Name, this_change(1))) = this_change(2);
        blupstats.Name(strcmp(blupstats.Name, this_change(1))) = this_change(2);
    end

    % prepare individual estimates for each parameter
    individual_estimates = NaN(nids, 1, npars);
    for pp = 1 : length(parnames2)
        tp = blupstats(strcmp(blupstats.Name, parnames(pp)), :);
        tp.Estimate = tp.Estimate + festats.Estimate(strcmp(festats.Name, parnames2(pp)));
        individual_estimates(:,:,pp) = tp.Estimate;
    end
    full_ind_estimates(:,tt,:) = individual_estimates;
    
    % export everything
    estimates(:,tt)   = festats.Estimate;
    lower_ci(:,tt)    = festats.Lower;
    upper_ci(:,tt)    = festats.Upper;
    t_stat(:,tt)      = festats.tStat;
    SE(:,tt)      = festats.SE;

    loglik(:,tt)      = rm.LogLikelihood;
    AIC(:,tt)         = rm.ModelCriterion{1,1};
    BIC(:,tt)         = rm.ModelCriterion{1,2};
    deviance(:,tt)    = rm.ModelCriterion{1,3};

    ordinary(:,tt)    = rm.Rsquared.Ordinary;
    adjusted(:,tt)    = rm.Rsquared.Adjusted;


end

rmdir(tmpout, 's')

%% Pack outputs
modelout.pars.estimates     = estimates;
modelout.pars.lower_ci      = lower_ci;
modelout.pars.upper_ci      = upper_ci;
modelout.pars.t_stat        = t_stat;
modelout.pars.parnames      = parnames;
modelout.infcriteria.loglik    = loglik;
modelout.infcriteria.AIC       = AIC;
modelout.infcriteria.BIC       = BIC;
modelout.infcriteria.deviance  = deviance;
modelout.rsqrd.ordinary     = ordinary;
modelout.rsqrd.adjusted     = adjusted;
modelout.full_ind_estimates = full_ind_estimates;
modelout.predictions.yhat   = predictions;
modelout.predictions.y      = averages;
modelout.residuals          = residuals;
modelout.stdz_residuals     = residuals ./ nanstd(residuals,0,1);
modelout.ids                = unique(data(:, idvar), 'stable');

end

% '(Intercept)'	'condition_congruent'	'condition_incongruent'	'rtwz'	'tonicbasez'	'trial'
function val = getOrDefault(s, field, defaultVal)
if isfield(s, field)
    val = s.(field);
else
    val = defaultVal;
end
end

%
% function averages = fit_and_export(data, tt, formula, tmpout,tsvar, yvar, idvar, verbose_fit)
%
% end