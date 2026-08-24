function aligned = tsglm_fit_one_point(tid, tsvar, formula, ref_names, npar, lkl, tt, ch)

tid_t = tid;

% parse the point to fit
if nargin < 8
    tid_t.(tsvar) = cellfun(@(x) x(tt), tid_t.(tsvar));
else
    tid_t.(tsvar) = cellfun(@(x) x(ch, tt), tid_t.(tsvar));
end

% fit the model
mdl = fitglm(tid_t, formula, 'Distribution', lkl);

% predictors names and betas
tmp_names = mdl.CoefficientNames;
tmp_betas = mdl.Coefficients.Estimate;

% align them to the prefixed order (fitglm may change the order of
% predictors in the output across different iterations)
aligned      = nan(npar, 1);
[lia, loc]   = ismember(ref_names, tmp_names);
aligned(lia) = tmp_betas(loc(lia));

% throw an error if you have any nan
if any(isnan(aligned)); error('nan estimates'); end


end