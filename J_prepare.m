format shortG

%
% STEP 1
%

% Load the SPM EEG file
D = spm_eeg_load('P10/spmeeg_P10BDF.mat');

% Load the events from CSV
event_data = readtable('P10/P10events_data.csv'); % Replace with the actual path

%
% STEP 2
%

% Jalynn:
% My understanding in trying to use only one condition at a time, is thatwe
% could only set the condition events for the condition we want to look at
% prior to even creating the trl file.

% Define the event values we need for each condition
condition_events = {
    [42, 12];  % Condition 1
    % [43, 13];  % Condition 2
    % [44, 14];  % Condition 3
    % [46, 16];  % Condition 4
};
condition_labels = { "Condition 2" };
% condition_labels = {"Condition 1", "Condition 2", "Condition 3", "Condition 4"};

%
% STEP 3
%

% Remove all channels except the first 3 and accel
channels_to_keep = [1, 2, 3];
accel_channels = find(contains(chanlabels(D), 'ACCEL'));
% Combine
channels_to_keep = [channels_to_keep, accel_channels];
% Create a new meeg object with only selected channels
S = struct();
S.D = D;
S.channels = D.chanlabels(channels_to_keep);
D_new = spm_eeg_crop(S);
D = D_new;

% Sampling frequency
fs = D.fsample;

%
% STEP 4
%

% Time window for epochs (in milliseconds)
prestim = 0; % Pre-stimulus time
poststim = 1000; % Post-stimulus time

% Initialize an empty trial structure
trials = [];
trl = [];
conditionlabels = [];

% Loop through each condition to define trial events
for cond = 1:length(condition_events)
    start_event = condition_events{cond}(1);
    end_event = condition_events{cond}(2);
    label = condition_labels{cond};
    
    % Filter the events for the current condition’s start markers
    cond_table = event_data(strcmp(event_data.type, 'Stimulus') & (event_data.value == start_event | event_data.value == end_event), :);

    if isempty(cond_table)
        continue
    end

    trials = [trials,
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', label),
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', append(label, "_TEST"))
    ];
    
    % Loop over each start event and define a trial
    for i = 1:2:height(cond_table)
        if i == height(cond_table)-1
            conditionlabels = [conditionlabels, append(label, "_TEST")];
        else
            conditionlabels = [conditionlabels, label];
        end
        trl = [trl; cond_table{i, "time"} cond_table{i+1, "time"} 0];
    end
end

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
% STEP 5
%

% Define the montage
% Keep only the channels specified earlier
montage_matrix = eye(D.nchannels);
% if D.nchannels = 4
% montage_matrix =
%     1     0     0     0
%     0     1     0     0
%     0     0     1     0
%     0     0     0     1
channels_to_remove = setdiff(1:D.nchannels, channels_to_keep);
% If D.nchannels = 4 and channels_to_keep = [1, 3], then
% channels_to_remove = setdiff([1, 2, 3, 4], [1, 3])
% channels_to_remove = [2, 4]
montage_matrix(channels_to_remove, :) = 0;
% montage_matrix =
%     1     0     0     0  % Channel 1 (retained)
%     0     0     0     0  % Channel 2 (removed)
%     0     0     1     0  % Channel 3 (retained)
%     0     0     0     0  % Channel 4 (removed)

% Apply the montage
S = struct();
S.D = D;
S.montage.tra = montage_matrix;
S.montage.labelorg = D.chanlabels;
S.montage.labelnew = D.chanlabels(channels_to_keep);
S.keepothers = 0;
D = spm_eeg_montage(S);

% Fix channel names
% Just naming them EEG1-3 for now since I don't know what they are suppose to
% be :)
for i = 1:3
    D = chanlabels(D, i, ['EEG' num2str(i)]);
end

% Rename Accelerometer channel?
for i = 1:length(accel_channels)
    accel_idx = accel_channels(i);
    D = chanlabels(D, accel_idx, ['ACC' num2str(i)]);
end

disp(trl);
disp(conditionlabels);
disp(trials);

trl = condition_timestamps;

%
% STEP 6 (Reduce TRL file)
% 

% remove all conditions except the 42/12
% My understanding is that this was a comment for an upcoming task - cut
% the trl file down to one condition only, starting with the LL condition
% I ended up modifying how we build trl file above instead of trying to
% remove or restart the trl file

%
% STEP 7 Jalynn: Converge - 1S Epoch
%

% Total recording duration in seconds
if ~isempty(condition_timestamps)
    recording_duration = condition_timestamps(:, 2) - condition_timestamps(:, 1);
    fprintf('Updated recording duration for Condition 42/12: %.2f seconds\n', recording_duration);
else
    recording_duration = D.nsamples / fs; % Default to total recording duration if condition timestamps are not set
end
% Define 1S epochs
% epoch_length = 1; % seconds! Could easily change to 2
% i'll be keeping the original trl file but will also be making an epoch 1S
% version of the trl file
% e_poststim = epoch_length * 1000;
% new trl for 1S epoch
e_trl = [];
% Loop through each row in the trl matrix
for i = 1:size(trl, 1)
    start_time = trl(i, 1);
    end_time = trl(i, 2);
    value = trl(i, 3);
    % Create 1-second intervals between start_time and end_time
    for t = start_time:1:end_time - 1  % Increment by 1 second, ensuring we don't exceed the end_time
        % Append the epoch to e_trl
        e_trl = [e_trl; t, t + 1, value];
    end
end
% Display a message if e_trl is still empty after processing
if isempty(e_trl)
    disp('No valid intervals were created. Check the trl data.');
else
    disp('New e_trl created successfully.');
end
% Save for reference
writematrix(e_trl, 'e_trl.csv');

if any(trl(:, 2) <= trl(:, 1))
    error('Invalid trl matrix: End times must be greater than start times.');
end
if any(trl(:) < 0) || any(trl(:, 2) > D.nsamples)
    error('Invalid trl matrix: Times must be within the data range.');
end

% Set up the SPM structure for epoching
S = struct();
S.D = D;
%S.trialdef = trials;
S.trl = trl;
S.conditionlabels = conditionlabels;
%S.timewin = [prestim, poststim]; % Define epoch time window in ms

% Create epochs based on the trials defined above
D = spm_eeg_epochs(S);

% Save the newly epoched data
D.save();

trl_1 = trl;
e_trl_1 = e_trl;
D_1 = D;