function aligned = tsglm_fit_one_point(tid, tsvar, formula, ref_names, npar, lkl, tt, ch)

% parse the point to fit
if nargin < 8
    tid.(tsvar) = cellfun(@(x) x(tt), tid.(tsvar));
else
    tid.(tsvar) = cellfun(@(x) x(ch, tt), tid.(tsvar));
end

% fit the model
mdl = fitglm(tid, formula, 'Distribution', lkl);

% predictors names and betas
tmp_names = mdl.CoefficientNames;
tmp_betas = mdl.Coefficients.Estimate;

% align them to the prefixed order
aligned      = nan(npar, 1);
[lia, loc]   = ismember(ref_names, tmp_names);
aligned(lia) = tmp_betas(loc(lia));

% throw an error if you have any nan
if any(isnan(aligned)); error('nan estimates'); end


end