%
% PREPARE.M
%

format shortG




%
% STEP 1: Load the EEG file & event data CSV
%

% Load the SPM EEG file
D = spm_eeg_load('C:/Users/Jalynn/Desktop/BDF_EEG_Files/P11/spmeeg_P11BDF.mat'); % MODIFY HERE (1/8)

% Load the events from CSV
event_data = readtable('C:/Users/Jalynn/Desktop/BDF_EEG_Files/P11/P11events_data.csv'); % MODIFY HERE (2/8)




%
% STEP 2: Fix channel names
%

% Fix channel names for the first 3 EEG channels
channel_names = {'Fz', 'Cz', 'Pz'};
for i = 1:3
    old_label = chanlabels(D, i); % Get the old channel label
    new_label = channel_names{i};
    D = chanlabels(D, i, new_label);
    disp(['Renamed channel ' old_label ' at index ' num2str(i) ' to ' new_label]);
end

% Find accelerometer channels
accel_channels = find(contains(chanlabels(D), 'ACC'));

% Rename Accelerometer channel

for i = 1:length(accel_channels)
    accel_idx = accel_channels(i);
    new_label = ['ACC' num2str(i)];
    D = chanlabels(D, accel_idx, new_label);
    disp(['Renamed accelerometer channel at index ' num2str(accel_idx) ' to ' new_label]);
end

% Find reference and ground channels
% TASK: 'REF' AND 'GND' ???
ref_channel = find(strcmp(chanlabels(D), 'REF')); % Adjust 'REF' to actual label if different
ground_channel = find(strcmp(chanlabels(D), 'GND')); % Adjust 'GND' to actual label if different

% Rename reference and ground channels
if ~isempty(ref_channel)
    D = chanlabels(D, ref_channel, 'REF');
    disp(['Renamed reference channel at index ' num2str(ref_channel) ' to REF']);
end

if ~isempty(ground_channel)
    D = chanlabels(D, ground_channel, 'GND');
    disp(['Renamed ground channel at index ' num2str(ground_channel) ' to GND']);
end




%
% STEP 3: Remove all channels except the first 3, accel, ref, and ground
%

% Identify the channels to keep
channels_to_keep = [1, 2, 3];

% Combine all channels to keep
channels_to_keep = unique([channels_to_keep, accel_channels, ref_channel, ground_channel]);

% Ensure indices are within the valid range
channels_to_keep = channels_to_keep(channels_to_keep <= D.nchannels);

% Create a new meeg object with only selected channels
S = struct();
S.D = D;
S.channels = D.chanlabels(channels_to_keep);
D = spm_eeg_crop(S);

% Update channels_to_keep to match the new dataset
channels_to_keep = 1:D.nchannels; % After cropping, indices reset to 1:N
disp(['Number of remaining channels: ' num2str(D.nchannels)]);
disp('Remaining channel labels:');
disp(D.chanlabels);





%
% STEP 4: Define the montage
%

% Define montage matrix
montage_matrix = eye(D.nchannels);

% Explanation:
% if D.nchannels = 4
% montage_matrix =
%     1     0     0     0
%     0     1     0     0
%     0     0     1     0
%     0     0     0     1

% Apply the montage to keep only the desired channels
channels_to_remove = setdiff(1:D.nchannels, channels_to_keep);

% Explanation:
% If D.nchannels = 4 and channels_to_keep = [1, 2, 3], then
% channels_to_remove = setdiff([1, 2, 3, 4], [1, 2, 3])
% channels_to_remove = [4]

% Actually removes channels
montage_matrix(channels_to_remove, :) = 0;

% Explanation:
% montage_matrix =
%     1     0     0     0  % Channel 1 (retained)
%     0     1     0     0  % Channel 2 (retained)
%     0     0     1     0  % Channel 3 (retained)
%     0     0     0     0  % Channel 4 (removed)

disp('Channels to keep:');
disp(channels_to_keep);

% Apply the montage
S = struct();
S.D = D;
S.montage.tra = montage_matrix;
S.montage.labelorg = D.chanlabels;
S.montage.labelnew = D.chanlabels(channels_to_keep);
S.keepothers = 0;
D = spm_eeg_montage(S);

disp('Applied montage successfully.');




%
% STEP 5: Define conditions for comparisons
%

% Define the event values we need for each condition
condition_events = { % MODIFY HERE (3/8)
    [42, 12];  % Condition 1 / LL
    % [43, 13];  % Condition 2 / LH
    % [44, 14];  % Condition 3 / HL
    % [46, 16];  % Condition 4 / HH
};

% Define the condition label
condition_labels = {"Condition 1"}; % MODIFY HERE (4/8)
% condition_labels = {"Condition 1", "Condition 2", "Condition 3", "Condition 4"};




%
% STEP 6: Fix the triggers
%

% Time window for epochs (in milliseconds)
prestim = 0; % Time before the stimulus
poststim = 1000; % Time after the stimulus

% Initialize an empty trial structure
trials = []; % Store trial definitions
trl = []; % Store trial timing (start, end, offset)
conditionlabels = string([]); % Initialize as an empty string array

% Loop through each condition to define trial events
for cond = 1:length(condition_events)
    start_event = condition_events{cond}(1);
    end_event = condition_events{cond}(2);
    label = string(condition_labels{cond});
    
    % Filter the events for the current condition’s start markers
    cond_table = event_data(strcmp(event_data.type, 'Stimulus') & (event_data.value == start_event | event_data.value == end_event), :);

    % Skip to the next condition if no relevant events are found
    if isempty(cond_table)
        continue
    end

    %
    % STEP 6.1: Define trial structure
    %

    if ~isempty(label)
        conditionlabels = [conditionlabels; label]; % Append if not empty
    else
        disp('Empty label detected, replacing with "Unknown!!!!".');
        conditionlabels = [conditionlabels; "Unknown"];
    end

    trials = [trials,
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', string([label])),
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', string([append(label, "_TEST")]))
    ];

    %
    % STEP 6.2: Generate Labels
    %
    
    % Loop over each start event and define a trial
    for i = 1:2:height(cond_table)
        if i == height(cond_table)-1
            conditionlabels = [conditionlabels, append(label, "_TEST")];
        else
            conditionlabels = [conditionlabels, label];
        end

        % Store trial start, end, and offset
        trl = [trl; cond_table{i, "time"} cond_table{i+1, "time"} 0];
    end
end

%
% STEP 6.3: Extract Timestamps
%

%% Get just the conditions, no trial level info
% Initialize a structure to store the results
condition_timestamps = [];

% Loop through each condition to extract timestamps
for cond = 1:length(condition_events)
    start_event = condition_events{cond}(1);
    end_event = condition_events{cond}(2);
    
    % Find the first occurrence of the start event and the last occurrence of the end event
    first_start_idx = find(event_data.value == start_event, 1, 'first');
    last_end_idx = find(event_data.value == end_event, 1, 'last');
    
    % Extract timestamps if both indices are valid
    if ~isempty(first_start_idx) && ~isempty(last_end_idx)
        first_start_time = event_data.time(first_start_idx);
        last_end_time = event_data.time(last_end_idx);
        
        % Append to condition_timestamps
        condition_timestamps = [condition_timestamps; first_start_time, last_end_time, 0];
    else
        fprintf('Events not found for Condition %d (%d-%d).\n', cond, start_event, end_event);
    end
end




%
% STEP 7: Epoch our data into 1 second intervals (to match the engagement paper)
%

fprintf('Starting epoching process...\n');

% Validate trl matrix
if any(trl(:, 2) <= trl(:, 1))
    error('Invalid trl matrix: End times must be greater than start times.');
end

if any(trl(:) < 0) || any(trl(:, 2) > D.nsamples)
    error('Invalid trl matrix: Times must be within the data range.');
end

% Sampling frequency
fs = D.fsample;

% Total recording duration in seconds
if ~isempty(condition_timestamps)
    recording_duration = condition_timestamps(:, 2) - condition_timestamps(:, 1);
    fprintf('Updated recording duration for Condition: %.2f seconds\n', recording_duration);
else
    error('No condition timestamps found. Cannot calculate recording duration.');
end

% new trl for 1 second epoching
epoch_trl = [];

% Generate 1 second epoch intervals
for i = 1:size(trl, 1)
    start_time = trl(i, 1);
    end_time = trl(i, 2);
    value = trl(i, 3);

    % Create 1-second intervals between start_time and end_time
    for t = start_time:1:end_time - 1  % Increment by 1 second, ensuring we don't exceed the end_time
        % Append the epoch to epoch_trl
        epoch_trl = [epoch_trl; t, t + 1, value];
    end
end

% Display a message if epoch_trl is still empty after processing
if isempty(epoch_trl)
    warning('No valid 1-second intervals were created. Verify trl data.');
else
    fprintf('Epoching completed successfully. Total epochs: %d\n', size(epoch_trl, 1));
end

% Save for reference
writematrix(epoch_trl, 'epoch_trl.csv');




%
% STEP 8: Baseline Correction
%

% Apply baseline correction to each epoch
% Define a structure to hold the baseline correction parameters
% S = struct();
S = [];
S.D = D;
S.timewin = [2762.35, 2777.488];  % Set baseline window (start and end samples)
D = spm_eeg_bc(S);  % Apply baseline correction

disp('Baseline correction applied successfully.');





%
% STEP 9: Output
%

disp(trl);
disp(conditionlabels);
disp(trials);


% Set up the SPM structure
S = struct();
S.D = D;
% This was here previously but I don't know if it's doing what we want it to
trl = condition_timestamps;
S.trl = trl;
S.conditionlabels = conditionlabels;
modified_D = spm_eeg_epochs(S);

modified_D.save();
modified_D_1 = modified_D; % MODIFY HERE (5/8)
trl_1 = trl; % MODIFY HERE (6/8)


% Set up the SPM structure for epoching
% epoch_S = struct();
% epoch_S.D = D;
% S.trl = epoch_trl;

% TASK: Jalynn added - warning about nonstring
% Doesn't help
% conditionlabels = [conditionlabels, string(label)];

% S.conditionlabels = conditionlabels;
% epoch_D = spm_eeg_epochs(S);

% epoch_D.save();
% epoch_D_1 = epoch_D; % MODIFY HERE (7/8)
% epoch_trl_1 = epoch_trl; % MODIFY HERE (8/8)
