



function tsglm_ssh_fit_eeg(run)
%% Paths 
proj_path = "/shared/projects/project_COBRA"; % SLURM path
% proj_path = "//129.199.81.50/cobra"; % local path

event = 'stimulus';
modname = event;

data_path = fullfile(proj_path, "Data", "DataPreproc");
ref_path  = fullfile(data_path, "eeg_glm", "ref_vars");
out_path  = fullfile(data_path, "eeg_glm", "estimates", modname);
eeg_path  = fullfile(data_path, "eeg_trim");

addpath(fullfile(proj_path, "Scripts", "Toolboxes", "tsglmm"));

% create a folder if it doesn't exist
if ~isfolder(out_path); mkdir(out_path); end

%% Model settings
n_runs              = 50; % number of runs (each run takes a certain number of particpants/channels)
expected_n_channels = 64; % number of channels

idvar   = "ppid";   % id Var
tsvar   = "eeg";    % time series var (eeg in this case)
lkl     = "Normal"; % likelihood ("Normal" for encoding)

% model formula
formula = "eeg ~ thr_cond * vol_cond * stochasticity_z + rt_wz + trial";

%% Load behavioral data 
data = readtable(fullfile(eeg_path, "behav_trim.csv"));
ids  = unique(data.ppid);

% normalise trial number and RT at the group level
data.trial = zscore(data.trial_num);
data.rt_z = zscore(data.rt_s);

% normalize RTs within participants
data.rt_wz = nan(height(data), 1);
for ii = 1 : length(ids)
    id_idx = find(strcmp(data.ppid, ids(ii)));
    data.rt_wz(id_idx) = zscore(data.rt_s(id_idx));
end

% make variables categorical
data.thr_cond = categorical(data.thr_cond);
data.thr_cond= reordercats(data.thr_cond, {'congruent','incongruent', 'neutral'});

data.vol_cond = categorical(data.vol_cond);
data.vol_cond= reordercats(data.vol_cond, {'stable','volatile'});

% define and normalise stochasticity
data.stochasticity_z = zscore(0.5-abs(data.probA-0.5));

% find EEG files...
eeg_files = dir(fullfile(eeg_path, ['*', event, '*.mat']));

% ... or throw an error if it fails
if isempty(eeg_files); error("No MAT files found in: %s", eeg_path); end

%% Load coefficient reference 
ref = load(fullfile(ref_path, "ref_names2.mat"));

if ~isfield(ref, "tmp_names")
    error("ref_names.mat does not contain a variable named ref_names.");
end

ref_names = ref.tmp_names;
ref_names = replace(ref_names, 'rt_z', 'rt_wz');
npar      = numel(ref_names);

%% Participant × channel combinations 
channels = 1:expected_n_channels;
[channel_grid, id_grid] = ndgrid(channels, 1:length(ids));
combos = [id_grid(:), channel_grid(:)];
n_combos = size(combos, 1);

% Divide combinations across SLURM runs. 
first_row = floor((run - 1) * n_combos / n_runs) + 1;
last_row  = floor(run * n_combos / n_runs);

rows_to_do = first_row:last_row;
nr = length(rows_to_do);

% display progress
fprintf( ...
    "Run %d: processing %d combinations, rows %d to %d.\n", ...
    run, nr, first_row, last_row ...
);

% kill the process if no rows to do
if nr == 0
    fprintf("Nothing assigned to run %d.\n", run);
    return
end

%% Fit the models assigned to this run
current_id = "";
for tr = 1 : nr

    % get the id and channel to fit
    row_number = rows_to_do(tr);

    idx = combos(row_number, 1); % the id
    id = ids(idx);

    ch  = combos(row_number, 2); % the channel
   
    % display progress
    fprintf( ...
        "   Run %d | combination %d/%d | participant %s | channel %d\n", ...
        run, tr, nr, char(id), ch ...
    );


    % load participant when participant changes otherwise keep the one
    % currently loaded
    if ~strcmp(id, current_id)
        
        % get ID
        current_id = id;
        tid = data(strcmp(data.(idvar), id), :);

        if isempty(tid)
            error("No behavioral rows found for participant %s.", id);
        end

        % find the related eeg data file
        t_eeg_file = eeg_files(cellfun(@(x) contains(x, current_id), {eeg_files.name})).name;
        
        % guard agaisnt multiple matches
        if numel(string(t_eeg_file)) ~= 1
            matching_names = string(t_eeg_file);
            error( ...
                "     Expected one EEG file for participant %s, but found %d: %s", ...
                id, n_matches, char(strjoin(matching_names, ", ")) ...
            );
        end
        
        % load eeg data
        loaded = load(fullfile(eeg_files(1).folder, t_eeg_file));

        % Accept either preprocessing variable name.
        if isfield(loaded, "data")
            eeg_data = loaded.data;
        else
            error("     File %s does not contain data field.", eeg_file.name);
        end
        
        % get activity by trial
        if ~isfield(eeg_data, "trial")
            error("     EEG structure in %s has no trial field.", eeg_file.name);
        end

        trials = eeg_data.trial;

        if ~iscell(trials)
            error("     eeg_data.trial must be a cell array.");
        end

        % ensure trials are nTrials × 1
        trials = trials(:);

        if height(tid) ~= numel(trials)
            error( ...
                ["     Trial-count mismatch for participant %s: " ...
                 "%d behavioral rows but %d EEG trials."], ...
                id, ...
                height(tid), ...
                numel(trials) ...
            );
        end

        first_trial_size = size(trials{1});

        if numel(first_trial_size) ~= 2
            error("     EEG trials must have dimensions channels × samples.");
        end

        nch   = first_trial_size(1); % number of channels
        tslen = first_trial_size(2); % length of the epochs

        same_trial_size = cellfun(@(x) isequal(size(x), first_trial_size), trials);

        if ~all(same_trial_size)
            error("     Not all EEG trials have the same dimensions for %s.", id);
        end

        if nch ~= expected_n_channels
            error( ...
                "     Expected %d channels for %s, but found %d.", ...
                expected_n_channels, id, nch ...
            );
        end

        tid.(tsvar) = trials;
    end


    % Fit every time sample ===============================================
    % Rows    = time samples
    % Columns = model coefficients
    pars = nan(tslen, npar);

    for tt = 1 : tslen

        % fit one sample point across all trials for this subject x channel
        % combo
        pars_one_sample = tsglm_fit_one_point( ...
            tid, ...
            tsvar, ...
            formula, ...
            ref_names, ...
            npar, ...
            lkl, ...
            tt, ...
            ch ...
        );

        % guard against non numeric outputs or outputs with wrong number of
        % elements
        if ~isnumeric(pars_one_sample) || numel(pars_one_sample) ~= npar
            error( ...
                ["     tsglm_fit_one_point returned an unexpected output " ...
                 "for participant %s, channel %d, sample %d. " ...
                 "Expected %d numeric coefficients."], ...
                id, ch, tt, npar ...
            );
        end
        
        % export the estimates for this sample point
        pars(tt, :) = pars_one_sample(:)';
    end


    % Save once, after fitting every sample
    if expected_n_channels < 100
        out_name = sprintf("%s_ch%.2d_%s.mat", char(id), ch, event);
    else
        out_name = sprintf("%s_ch%.3d_%s.mat", char(id), ch, event);
    end

    save( ...
        fullfile(out_path, out_name), ...
        "pars", ...
        "ref_names", ...
        "id", ...
        "ch" ...
    );

end


end