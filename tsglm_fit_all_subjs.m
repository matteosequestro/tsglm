function [modelout] = tsglm_fit_all_subjs(data, formula, cfg)

% parse variables from cfg
idvar = cfg.idvar;


lkl = cfg.glm_likelihood;

% parse ids
ids   = unique(data{:,idvar});
nids  = numel(ids);

% parse time-series variable's name
tsvar = cfg.tsvar;

if ~isfield(cfg, 'eeg_path') % if it's not eeg
 
    % get time series variable and it size
    size_ts = size(data.(tsvar){1});

    % length of the time series
    tslen = size_ts(2);

    % build global reference model (the order of predictors from fitglm may
    % change across runs, which will be a problem, so I fit a first one and
    % force this order to the following)
    data_tmp = data;
    data_tmp.(tsvar) = cellfun(@(x) x(1), data_tmp.(tsvar));

    ref_fit   = fitglm(data_tmp, formula, 'Distribution', lkl);
    ref_names = ref_fit.CoefficientNames;    % coefficients names (i.e., the order to follow)
    npar      = numel(ref_names);            % number of coefficients

    % preallocate estimates array
    pars_series = nan(nids, tslen, npar);

    % main loop
    for ii = 1:nids
        % parse participant
        tid = data(strcmp(data{:,idvar}, ids{ii}),:);

        % loop through samples
        parfor tt = 1 : tslen
            % fit export to the preallocated array
                pars_series(ii, tt, :) = tsglm_fit_one_point(tid, tsvar, formula, ref_names,npar, lkl, tt);
        end
    end

    
else % if it's eeg
    % get eeg data files
    eeg_files = dir(fullfile(cfg.eeg_path, '*.mat'));

    % create a fake dependent variable to run the mock model
    data_tmp = data;
    data_tmp.(tsvar) = normrnd(0, 1, height(data_tmp), 1);

    % fit the mock model
    ref_fit   = fitglm(data_tmp, formula, 'Distribution', lkl);
    ref_names = ref_fit.CoefficientNames;    % coefficients names (i.e., the order to follow)
    npar      = numel(ref_names);            % number of coefficients


    % main loop
    for ii = 1 : nids
        % parse participant
        id = ids{ii};
        tid = data(strcmp(data{:,idvar}, id),:);

        % get eeg data for the subject
        which_eeg = find(cellfun(@(x) contains(x, id), {eeg_files.name}));
        ts0 = load(fullfile(eeg_files(which_eeg).folder, eeg_files(which_eeg).name));
        tid.(tsvar) = ts0.data.trial';

        % length of the time series
        size_ts = size(tid.(tsvar){1});
        tslen = size_ts(2);
        nch = size_ts(1); % usually > 1 only with eeg

        % preallocate the parameters array
        if ~exist("pars_series", "var")
            pars_series = nan(nids, nch, tslen, npar);
        end

        % loop through channels ----- 
        for ch = 1 : nch
            % show progress if necessary
            if cfg.verbose; fprintf('subj: %d/%d, nids | channel: %d/%d\n', ii, nids, ch, nch); end
            
            % loop through samples
            parfor tt = 1 : tslen
                % fit export to the preallocated array
                pars_series(ii, ch, tt, :) = tsglm_fit_one_point(tid, tsvar, formula, ref_names,npar, lkl, tt, ch);
            end % ----- sample loop
        end % ----- channel loop
    end  % ----- subj loop

end



% output
modelout.pars_series = pars_series;
modelout.par_names   = ref_names;
modelout.ids         = ids;

end