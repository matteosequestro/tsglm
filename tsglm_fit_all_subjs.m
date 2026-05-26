function [modelout] = tsglm_fit_all_subjs(data, formula, glm_likelihood, cfg)

% parse variables from cfg
idvar = cfg.idvar;
tsvar = cfg.tsvar;

% parse ids
ids   = unique(data{:,idvar});
nids  = numel(ids);

% length of the time series
tslen = numel(data.(tsvar){1});

% build global reference model (the order of predictors from fitglm may
% change across runs, which will be a problem, so I fit a first one and
% force this order to the following)
data_tmp = data;
data_tmp.(tsvar) = cellfun(@(x) x(1), data_tmp.(tsvar));

ref_fit   = fitglm(data_tmp, formula, 'Distribution', glm_likelihood);
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

        % parse the sample to fit
        tad = tid;
        tad.(tsvar) = cellfun(@(x) x(tt), tad.(tsvar));

        % fit the model
        mdl = fitglm(tad, formula, 'Distribution', glm_likelihood);
        
        % predictors names and betas
        tmp_names = mdl.CoefficientNames;
        tmp_betas = mdl.Coefficients.Estimate;
        
        % align them to the prefixed order
        aligned      = nan(npar,1);
        [lia, loc]   = ismember(ref_names, tmp_names);
        aligned(lia) = tmp_betas(loc(lia));

        % throw an error if you have any nan
        if any(isnan(aligned)); error('nan estimates'); end

        % export to the preallocated array
        pars_series(ii, tt, :) = aligned;
    end
end



% output
modelout.pars_series = pars_series;
modelout.par_names   = ref_names;
modelout.ids         = ids;

end