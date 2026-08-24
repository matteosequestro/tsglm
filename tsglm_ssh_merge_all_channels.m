clear all; clc; close all

%% Paths 
% proj_path = "/shared/projects/project_COBRA";
proj_path = "//129.199.81.50/cobra";

data_path   = fullfile(proj_path, "Data", "DataPreproc");
eeg_path    = fullfile(data_path, "eeg_trim");
glm_in_path = fullfile(data_path, "eeg_glm", "estimates", "stimulus");
out_path    = fullfile(data_path, "eeg_glm", "ready");

if ~isfolder(out_path)
    mkdir(out_path);
end

%% Input files 
glm_files = dir(fullfile(glm_in_path, "*ch*_stimulus.mat"));
eeg_files = dir(fullfile(eeg_path, "*stimulus_trim.mat"));

if isempty(glm_files)
    error("No GLM files found in %s.", glm_in_path);
end

if isempty(eeg_files)
    error("No EEG files found in %s.", eeg_path);
end

modname = "stimulus";

%% Derive participant IDs from filenames such as ID_ch-01.mat 
glm_names = string({glm_files.name});

if any(~contains(glm_names, "_ch"))
    error("Some GLM filenames do not contain the expected '_ch-' token.");
end

ids = unique(extractBefore(glm_names, "_ch"), "stable");
nn  = numel(ids);

%% Read coefficient information from the first GLM file
first_glm = load(fullfile(glm_files(1).folder, glm_files(1).name));

if ~isfield(first_glm, "ref_names") || ~isfield(first_glm, "pars")
    error("The first GLM file does not contain ref_names and pars.");
end

original_par_names = string(first_glm.ref_names(:));
n_pars             = numel(original_par_names);

% Produce valid and unique MATLAB structure-field names
par_names = matlab.lang.makeValidName(original_par_names);
par_names = matlab.lang.makeUniqueStrings(par_names);

% Preserve the mapping between original coefficient names and field names
coefficient_info = table( ...
    original_par_names, ...
    string(par_names(:)), ...
    'VariableNames', {'original_name', 'field_name'} ...
);

%% Preallocate output structure 
all_coefs = struct();

for pp = 1:n_pars
    all_coefs.(par_names{pp}) = cell(nn, 1);
end

reference_labels = [];
reference_time   = [];

%% Loop through participants
for ii = 1 : nn
 
    id = ids(ii);

    fprintf("%d/%d | participant %s\n", ii, nn, id);

    % Find this participant's GLM files using an exact filename prefix
    glm_match = startsWith(glm_names, id + "_ch");
    glm_i     = glm_files(glm_match);

    if isempty(glm_i)
        error("No GLM files found for participant %s.", id);
    end

    % Find the participant's EEG file as a complete filename token
    eeg_t = eeg_files( ...
        cellfun(@(x) contains(x, id), {eeg_files.name}) ...
        ).name;
    
    % Ensure that the number of files is correct
    if isempty(eeg_t);            error("no EEG found for %s", id);          end
    if length(string(eeg_t)) > 1; error("multiple matching EEG for %s", id); end
    
    % Load and validate EEG data
    loaded   = load(fullfile(eeg_files(1).folder, eeg_t));

    if isfield(loaded, "data")
        eeg_data = loaded.data;
    else
        error("EEG file %s does not contain data.", eeg_t);
    end

    % Validate EEG metadata
    if ~isfield(eeg_data, "label")
        error("EEG file %s has no label field.", eeg_file.name);
    end

    if ~isfield(eeg_data, "time") || isempty(eeg_data.time)
        error("EEG file %s has no valid time field.", eeg_file.name);
    end
    
    % Validate channel info
    labels     = eeg_data.label(:);
    n_channels = numel(labels);

    if iscell(eeg_data.time)
        time_vector = eeg_data.time{1};
        same_time = cellfun(@(x) isequal(x, time_vector), eeg_data.time);

        if ~all(same_time)
            error("Not all trials have the same time axis for participant %s.", id);
        end
    else
        time_vector = eeg_data.time;
    end

    time_vector = time_vector(:)';
    n_time      = numel(time_vector);

    % Ensure that there is a GLM estimate file for every channel
    if numel(glm_i) ~= n_channels
        error( ...
            "Participant %s has %d channel GLM files but %d EEG channel labels.", ...
            id, numel(glm_i), n_channels ...
        );
    end

    % Ensure channel order and time axis agree across participants
    if ii == 1
        reference_labels = labels;
        reference_time   = time_vector;
    else
        if ~isequal(labels, reference_labels)
            error("Channel labels/order differ for participant %s.", id);
        end

        if ~isequal(time_vector, reference_time)
            error("Time axis differs for participant %s.", id);
        end
    end

    % Construct a clean FieldTrip timelock template
    template = struct();

    template.label              = labels;
    template.time               = time_vector;
    template.dimord             = "chan_time";

    template.cfg                = struct();
    template.cfg.source_file    = eeg_t;

    % Ensure that the template contains an elec field
    if isfield(eeg_data, "elec")
        template.elec = eeg_data.elec;
    else
        error("eeg_data should contain an 'elec' field but it doesn't for participant %s", id)
    end

    % Preallocate one field per coefficient
    for pp = 1:n_pars
        all_coefs.(par_names{pp}){ii} = template;
        all_coefs.(par_names{pp}){ii}.avg = nan(n_channels, n_time);
    end

    seen_channels = false(n_channels, 1);

    % Loop through glm files to extract estimates
    for gg = 1 : numel(glm_i)
        
        % load glm
        glm_t = load(fullfile(glm_i(gg).folder, glm_i(gg).name));

        % ensure that the required fields are present 
        required_fields = ["pars", "ref_names", "ch"];

        for ff = 1:numel(required_fields)
            if ~isfield(glm_t, required_fields(ff))
                error( ...
                    "GLM file %s lacks field %s.", ...
                    glm_i(gg).name, required_fields(ff) ...
                );
            end
        end
        
        % get and validate current channel
        ch = glm_t.ch;

        if ~isscalar(ch) || ch < 1 || ch > n_channels 
            error( ...
                "Invalid channel index in GLM file %s: %s.", ...
                glm_i(gg).name, mat2str(ch) ...
            );
        end

        if seen_channels(ch)
            error( ...
                "Duplicate estimate for participant %s, channel %d.", ...
                id, ch ...
            );
        end

        seen_channels(ch) = true;

        % validate length of the time series of coefficients
        if size(glm_t.pars, 1) ~= n_time
            error( ...
                ["Time dimension mismatch in %s: expected %d samples, " ...
                 "found %d."], ...
                glm_i(gg).name, n_time, size(glm_t.pars, 1) ...
            );
        end
    
        % validate the number of coefficients/parameters
        if size(glm_t.pars, 2) ~= n_pars
            error( ...
                ["Coefficient-count mismatch in %s: expected %d, " ...
                 "found %d."], ...
                glm_i(gg).name, n_pars, size(glm_t.pars, 2) ...
            );
        end
    
        % validate coefficients names
        this_ref_names = string(glm_t.ref_names(:));
        if ~isequal(this_ref_names, original_par_names)
            error( ...
                "Coefficient names/order differ in GLM file %s.", ...
                glm_i(gg).name ...
            );
        end
        
        % put coefficients in the preallocated array
        for pp = 1:n_pars
            all_coefs.(par_names{pp}){ii}.avg(ch, :) = glm_t.pars(:, pp)';
        end
    end

    if ~all(seen_channels)
        missing_channels = find(~seen_channels);

        error( ...
            "Participant %s is missing channels: %s", ...
            id, mat2str(missing_channels') ...
        );
    end
end

%% Export
out_name = modname + ".mat";

save( ...
    fullfile(out_path, out_name), ...
    "all_coefs", ...
    "coefficient_info", ...
    "ids", ...
    "-v7.3" ...
);