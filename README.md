# TSGLM

**Time-series Generalized Linear Models with cluster-based permutation testing in MATLAB**

TSGLM is a MATLAB <u>*(work in progress)*</u> toolbox for fitting generalized linear models (GLMs) to trial-wise time-series data and performing group-level cluster-based permutation tests on the resulting coefficient time courses.

It was developed for the analysis of psychophysiological and neurophysiological signals in experimental data, where the effect of trial-level predictors is estimated independently at each time point and statistical inference must account for multiple comparisons across time.

Typical applications include time-resolved analyses of:

* EEG-derived signals
* pupil responses
* ECG / heart-rate responses
* other continuous physiological or behavioural time series

The toolbox combines MATLAB's `fitglm` framework with FieldTrip's cluster-based permutation procedures.

## Overview

The basic analysis pipeline is:

1. For each participant, fit the same GLM independently at every time point of the signal.
2. Store the resulting regression coefficients as participant × time × parameter arrays.
3. Test each coefficient time series against zero at the group level.
4. Identify temporally contiguous clusters of effects using cluster-based permutation testing.
5. Extract participant-level mean coefficients within significant clusters for further analyses or visualization.

This allows conventional trial-level regression models to be applied to continuous time-series responses while controlling for multiple comparisons across time.

## Installation

Clone or download this repository and add it, including the `utils` directory, to the MATLAB path:

```matlab
addpath('/path/to/tsglm')
addpath('/path/to/tsglm/utils')
```

FieldTrip must also be installed and available on the MATLAB path.

## Data format

Input data should be provided as a MATLAB table with **one row per trial**.

The table must contain:

* a participant identifier;
* one column for each predictor included in the GLM;
* one cell-array column containing the time series for each trial.

For example:

```matlab
data =

    id      condition     prediction_error      signal
    __      _________     ________________      __________

    "01"       0               -0.42             {1x500 double}
    "01"       1                0.81             {1x500 double}
    "01"       0                0.15             {1x500 double}
    "02"       1               -0.30             {1x500 double}
    ...
```

All trial-level time series must have the same number of samples.

## Basic usage

Define the relevant columns:

```matlab
cfg.idvar = 'id';
cfg.tsvar = 'signal';
```

Specify a model using MATLAB/Wilkinson formula notation:

```matlab
formula = 'signal ~ 1 + prediction_error * condition';
```

Then run:

```matlab
[modelout, clust_par_means] = ...
    tsglm_run_clustperm(data, formula, cfg);
```

The toolbox will:

* fit the GLM separately for every participant and time point;
* construct coefficient time series;
* perform cluster-based permutation tests for each coefficient;
* return information about detected clusters;
* optionally plot the group-level coefficient time courses.

## Configuration

The main configuration structure is `cfg`.

### Required fields

```matlab
cfg.idvar
```

Name of the participant-ID column.

```matlab
cfg.tsvar
```

Name of the table column containing the trial-level time series.

### GLM likelihood

The default model is Gaussian:

```matlab
cfg.glm_likelihood = 'Gaussian';
```

Other likelihoods supported by MATLAB's `fitglm` can be specified when appropriate.

### FieldTrip permutation settings

FieldTrip options can be supplied through:

```matlab
cfg.fieldtrip_cfg
```

For example:

```matlab
cfg.fieldtrip_cfg.numrandomization = 5000;
cfg.fieldtrip_cfg.alpha = 0.05;
cfg.fieldtrip_cfg.clusteralpha = 0.05;
```

Default settings include:

```matlab
method           = 'montecarlo'
statistic        = 'depsamplesT'
correctm         = 'cluster'
clusterstatistic = 'maxsum'
tail             = 0
clustertail      = 0
clusteralpha     = 0.05
alpha            = 0.05
```

The analysis is one-dimensional in time, so no spatial-neighbour structure is required.

### Plotting

Coefficient time courses and detected clusters are plotted by default.

Plotting can be disabled with:

```matlab
cfg.wantplot_perm = 0;
```

A stimulus/event onset can also be marked in the time-series plots:

```matlab
cfg.stim_onset_time = 100;
```

where the value corresponds to the relevant sample index.

## Outputs

### `modelout`

The main output structure contains the participant-level coefficient time series and results of the cluster analysis.

Important fields include:

```matlab
modelout.pars_series
```

A three-dimensional array:

```text
participants × time points × regression coefficients
```

```matlab
modelout.par_names
```

Names of the GLM coefficients.

```matlab
modelout.ids
```

Participant identifiers.

```matlab
modelout.obs_clusters_sum
```

Summary of clusters detected for the different regression coefficients, including their temporal extent and permutation statistics.

### `clust_par_means`

A table containing, for each participant, the mean regression coefficient within each detected cluster.

This output can be used for visualization, correlations with individual-difference variables, or other secondary analyses.

## Main functions

### `tsglm_run_clustperm`

Main analysis function. Fits the time-resolved GLMs and performs cluster-based permutation testing.

```matlab
[modelout, clust_par_means, outputText, outputStats] = ...
    tsglm_run_clustperm(data, formula, cfg);
```

A previously fitted `modelout` can also be supplied to rerun the permutation stage without refitting all participant-level GLMs.

### `tsglm_fit_all_subjs`

Fits a GLM separately at every time point for every participant.

```matlab
modelout = tsglm_fit_all_subjs(data, formula, glm_likelihood, cfg);
```

### `tsglm_clusters_parmeans`

Extracts each participant's mean coefficient within detected clusters:

```matlab
clust_par_means = tsglm_clusters_parmeans(modelout);
```

### `tsglm_plot_estimates`

Plots group-average coefficient time courses with confidence intervals and detected clusters:

```matlab
tsglm_plot_estimates(modelout, cfg);
```

### `tsglm_plot_individual_pars`

Plots the coefficient time courses of individual participants:

```matlab
tsglm_plot_individual_pars(modelout);
```

### `tsglm_plot_clust_par_means`

Visualizes participant-level cluster averages and computes one-sample Bayes factors for cluster-level effects.

### `tsglm_write_stats`

Produces a compact summary table of detected clusters and their statistics.

## Example

A typical analysis might look like:

```matlab
% Add toolbox
addpath('/path/to/tsglm')
addpath('/path/to/tsglm/utils')

% Configure data
cfg.idvar = 'subject';
cfg.tsvar = 'pupil';

% Number of permutations
cfg.fieldtrip_cfg.numrandomization = 5000;

% Statistical thresholds
cfg.fieldtrip_cfg.clusteralpha = 0.05;
cfg.fieldtrip_cfg.alpha = 0.05;

% Example trial-level model
formula = 'pupil ~ 1 + prediction_error * threat + trial + reaction_time';

% Run analysis
[modelout, cluster_means] = ...
    tsglm_run_clustperm(data, formula, cfg);

% Display cluster statistics
stats = tsglm_write_stats(modelout);

% Plot individual coefficient trajectories
tsglm_plot_individual_pars(modelout);
```

For each participant and each time point, the toolbox fits the specified model to the trials:

```text
signal(t) ~ prediction_error × threat + trial + reaction_time
```

This produces a separate time course for the intercept, each main effect, and each interaction. Cluster-based permutation testing is then performed independently on these coefficient time courses.

## Additional utilities

The repository also contains supplementary plotting and analysis functions in:

```text
extra_plotting/
```

including utilities for correlation plots and within-subject error estimates.

General helper functions used by the main analysis pipeline are stored in:

```text
utils/
```

## Dependencies

Required:

* MATLAB
* Statistics and Machine Learning Toolbox
* Parallel Computing Toolbox for parallelized fitting
* FieldTrip

Some optional plotting/statistical functions may require additional MATLAB packages or functions.

## Statistical approach

Cluster-based permutation testing is based on the approach described by:

> Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. *Journal of Neuroscience Methods, 164*(1), 177–190.

The procedure controls the multiple-comparison problem arising from testing effects at many consecutive time points while retaining sensitivity to temporally extended effects.

Importantly, cluster-level inference supports conclusions about the presence of an effect across a temporal cluster, rather than precise inference about the onset or duration of that effect.


## License

This repository is currently distributed under the Creative Commons Attribution-NonCommercial 4.0 International license (CC BY-NC 4.0).
