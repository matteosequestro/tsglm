

function [modelout, clust_par_means, outputText, outputStats] = tsglm_run_clustperm(data, formula, cfg, modelout)
% ----------------------------------------------------------------------------------------------------------------------------------------------
% Performs cluster based correction on parameter estimates from 
% generalised linear model on a time-series
%
% INPUTS:
% data                  : full dataset having at least a column for ID,
%                           one column for each predictor in your model,
%                           and one column with the time series in a cell.
%                           Each row is a trial.
% formula               : the model formula in Wilkinson notation (i.e., R-like
%                           notation, e.g., 'y ~ 1 + x1 * x2 + (1 + x1 * x2 | id)'.
%                           It must be a string.
% cfg.niter             : number of permutations.
% cfg.wantplot_perm     : 1 if you want to plot time-series for each
%                           parameter with significant (blue) and
%                           non-significant (gray) clusters as shaded
%                           reptangles. 0 if you don't want such a plot.
% cfg.perm_alpha        : permutation alpha. Montecarlo p-values below such
%                           a value will be considered as significant.
% cfg.want_parallel_fit : 1 will fi the models across different time-points
%                           in parallel (much faster but computationally
%                           more expensive); 0 will do one point at the
%                           time.
% OUTPUTS:
% modelout              : ouput of the model on true dataset (i.e., not
%                           permuted. You will get a "pars" field
%                           (parameter estimate, lower and upper confidence
%                           interval, tvalues), a "criteria" field with
%                           log-likelihood of the model, as well as AIC,
%                           BIC and deviance criteria (for each sample
%                           point; and an "rsqrd" field with ordinary andf$costsd
%                           adjusted R-squared
% obs_clusters_sum      : output of the cluster based permutation. For each
%                           coefficient in the model and each cluster
%                           detected on the true dataset (i.e., not
%                           permutedyou), get a monte-carlo
%                           p-value, first sample in the cluster, length
%                           of the cluster and t-mass.
% ----------------------------------------------------------------------------------------------------------------------------------------------


% Set defaults if not set in the cfg
glm_likelihood                  = get_or_default(cfg, 'glm_likelihood', 'Gaussian');
cut_short_clusters              = get_or_default(cfg, 'cut_short_clusters', 0);
minlength                       = get_or_default(cfg, 'minlength', 0);
wantplot_perm                   = get_or_default(cfg, 'wantplot_perm', 1);
modname                         = get_or_default(cfg, 'modname', "");               % just for plotting purposes, e.g., if you need model comparison
fieldtrip_cfg                   = get_or_default(cfg, 'fieldtrip_cfg', struct()); 

% Take time
tic

% Specific fields for the cluster-based permutation, probably you won't
% need to change them.
fieldtrip_cfg.uvar       = 1;  % subject
fieldtrip_cfg.ivar       = 2;  % condition
fieldtrip_cfg.spmversion = 'spm12';
fieldtrip_cfg.feedback   = 'no';

% Get defaults for fieldtrip
fieldtrip_cfg.numrandomization  = get_or_default(fieldtrip_cfg, 'numrandomization', 1e3);              % number of permutations
fieldtrip_cfg.method            = get_or_default(fieldtrip_cfg, 'method',           'montecarlo');     % use monte carlo sampling for 'p-values'
fieldtrip_cfg.statistic         = get_or_default(fieldtrip_cfg, 'statistic',        'depsamplesT');    % paired t-test (equivalent to t-test agains 0 if contrasted with a null matrix of zeros)
fieldtrip_cfg.correctm          = get_or_default(fieldtrip_cfg, 'correctm',         'cluster');        % correction method
fieldtrip_cfg.clusterstatistic  = get_or_default(fieldtrip_cfg, 'clusterstatistic', 'maxsum');         % statistic to compute on permutations
fieldtrip_cfg.tail              = get_or_default(fieldtrip_cfg, 'tail',             0);                % two-sided test to find clusters
fieldtrip_cfg.clustertail       = get_or_default(fieldtrip_cfg, 'clustertail',      0);                % compare observed clusters against positive and negative permutation clusters
fieldtrip_cfg.clusteralpha      = get_or_default(fieldtrip_cfg, 'clusteralpha',     0.05);             % alpha for cluster detection
fieldtrip_cfg.alpha             = get_or_default(fieldtrip_cfg, 'alpha',            0.05);             % alpha for monte carlo testing
fieldtrip_cfg.neighbours        = get_or_default(fieldtrip_cfg, 'neighbours',       []);               % no spatial neighbours since it's unidimensional
fieldtrip_cfg.correcttail       = get_or_default(fieldtrip_cfg, 'correcttail',      'alpha');


% Fit the model
if nargin < 4
    disp('_________________________________________________')
    disp(['fitting a ', change_text(glm_likelihood), ' model with the formula: ', formula '...'])
   
    % fit the model on everyone
    tic
    modelout = tsglm_fit_all_subjs(data, formula, glm_likelihood, cfg);
    fprintf('fitting done in %.2fs!\n', round(toc,2));
    
    % Extract some variables from the model and add them modelout
    modelout.formula    = formula;
    modelout.modname    = modname;
    modelout.likelihood = glm_likelihood; 
end

% add some info to modelout
npar     = size(modelout.pars_series, 3); % number of parameters (a.k.a., coefficients)
tslen    = size(modelout.pars_series, 2); % length of the time-series
parnames = modelout.par_names;

% Unpack individual parameter estimates
pars_series = modelout.pars_series;

% Run Fieldtrip cluster based permutationn
n_subjs = size(pars_series, 1);
timevec = 1 : size(pars_series, 2);

% Design matrix for fieldtrip
design(1,:) = repmat(1 : n_subjs, 1, 2);
design(2,:) = [ones(1, n_subjs), ones(1, n_subjs) .* 2];
fieldtrip_cfg.design = design;

% Preallocate output arrays
outputText = cell(npar, 1);
outputStats = cell(npar, 1);

modelout.alpha = fieldtrip_cfg.alpha;

% Loop through parameters
disp('permuting...')
for pp = 1 : size(pars_series, 3)
    mySignal = pars_series(:,:,pp);

    % Shape parameters in fieldtrip format
    subjdata = cell(1, n_subjs);
    for subj = 1 : n_subjs
        subjdata{subj}.avg    = mySignal(subj,:);   % 1 x nTime of signal
        subjdata{subj}.time   = timevec;            % 1 x nTime 1 to last moment
        subjdata{subj}.label  = {'signal'};         % single channel label
        subjdata{subj}.dimord = 'chan_time';
    end

    % Create contrast data (only once, then it's the same for all pars)
    if pp == 1
        zero_data = subjdata;
        for dd = 1 : length(zero_data)
            zero_data{dd}.avg = zero_data{dd}.avg .* 0;
        end
    end

    % Run cluster based permutation (Maris and Oostenveld, 2007, Journal of Neuroscience Methods)
    outputText{pp} = evalc('stat = ft_timelockstatistics(fieldtrip_cfg, subjdata{:}, zero_data{:});');
    outputStats{pp} = stat;

    % Extract positive cluster info
    if isfield(stat, "posclusters")
        n_pclusters = length(stat.posclusters);
        for cl = 1 : n_pclusters
            clus_pos = find(stat.posclusterslabelmat == cl);
            stat.posclusters(cl).first      = clus_pos(1);
            stat.posclusters(cl).last       = clus_pos(end);
            stat.posclusters(cl).length     = length(clus_pos);
        end
        posT = struct2table(stat.posclusters);
    else
        posT = [];
    end

    % Extract negative cluster info
    if isfield(stat, "negclusters")
        n_nclusters = length(stat.negclusters);
        for cl = 1 : n_nclusters
            clus_neg                        = find(stat.negclusterslabelmat == cl);
            stat.negclusters(cl).first      = clus_neg(1);
            stat.negclusters(cl).last       = clus_neg(end);
            stat.negclusters(cl).length     = length(clus_neg);
        end
        negT = struct2table(stat.negclusters);
    else
        negT = [];
    end


    % Merge positive and negative clusters
    if ~isempty(posT) && ~isempty(negT)
        clusters = [posT; negT];
    elseif ~isempty(posT)
        clusters = posT;
    elseif ~isempty(negT)
        clusters = negT;
    else
        continue  % skips to next iteration of the for-loop
    end
    clear posT negT

    clusters            = sortrows(clusters, 'first');                      % Sort clusters by timing
    clusters.clustern   = (1 : height(clusters))';                          % Add cluster number
    clusters.pred       = repmat(change_text(parnames{pp}), height(clusters), 1);        % Add predictor name

    % Export to the main dataset (or create it)
    if ~exist("obs_clusters_sum", "var")
        obs_clusters_sum = clusters;
    else
        obs_clusters_sum = [obs_clusters_sum; clusters];
    end
end



% Cut short clusters
if cut_short_clusters
    obs_clusters_sum(obs_clusters_sum.length < minlength,:) = [];
end

% Export observed cluster summary
modelout.obs_clusters_sum = obs_clusters_sum;

% Plots
if wantplot_perm
    tsglm_plot_estimates(modelout, cfg);
end 

% Compute mean parameter estimates for each individual and significant cluster
if nargout > 1 % Do this only if it's required by the output to save time
    clust_par_means = tsglm_clusters_parmeans(modelout); % the function could have a better name
end

% Write stats if required
if cfg.write_stats
    tsglm_write_stats(modelout);
end

% Display total elapsed time
fprintf('done in %.2fs!\n', round(toc,2));


end
