%
% 12.19.2024
% I just had the biggest aha moment. We only want one line in the trl - the
% start and end of one condition like 42/12. So why did we want 7 labels
% in the conditionlabels? This created an unequal amount of information in
% the S variable creating a line of issues starting with string being set
% to unknown, going down a path and eventually coming back to the same
% conclusion that if we wanted to use all 7 labels, we would have to use
% all 7 sets of cycles creating 7 unequal sets of trials and spm won't work
% with datasets that are unequal hence the previous error: all trials should
% have identical and positive lengths. 

% What does this mean for us now. The outputs are all correct right now
% and create a one line trl file with one label using conditional_label
% from the beginning instead of the heavily modified conditionlabels. This
% means that a lot of the code in this script is unnecessary and should
% probably be cut. This is now an upcoming task as the goal moving forward
% is just to have working code. 

% Also note that it sounds like using fieldtrip would have allowed for
% variable-length epochs

%
% PREPARE.M
%

format shortG




%
%% STEP 1: Load the EEG file & event data CSV
%

% Load the SPM EEG file
D_original = spm_eeg_load('C:/Users/Jalynn/Desktop/BDF_EEG_Files/P11/spmeeg_P11BDF.mat'); % MODIFY HERE (1/13)

% Load the events from CSV
event_data = readtable('C:/Users/Jalynn/Desktop/BDF_EEG_Files/P11/P11events_data.csv'); % MODIFY HERE (2/13)




%
%% STEP 2: Fix channel names
%

% Fix channel names for the first 3 EEG channels
channel_names = {'Fz', 'Cz', 'Pz'};
for i = 1:3
    old_label = chanlabels(D_original, i); % Get the old channel label
    new_label = channel_names{i};
    D_original = chanlabels(D_original, i, new_label);
    disp(['Renamed channel ' old_label ' at index ' num2str(i) ' to ' new_label]);
end

% Find accelerometer channels
accel_channels = find(contains(chanlabels(D_original), 'ACC'));

% Rename Accelerometer channel
for i = 1:length(accel_channels)
    accel_idx = accel_channels(i);
    new_label = ['ACC' num2str(i)];
    D_original = chanlabels(D_original, accel_idx, new_label);
    disp(['Renamed accelerometer channel at index ' num2str(accel_idx) ' to ' new_label]);
end

% Find reference and ground channels
% TASK: 'REF' AND 'GND' ???
% ref_channel = find(strcmp(chanlabels(D_original), 'REF')); % Adjust 'REF' to actual label if different
ref_channel = find(strcmp(chanlabels(D_original), 'F7')); % Adjust 'REF' to actual label if different

% ground_channel = find(strcmp(chanlabels(D_original), 'GND')); % Adjust 'GND' to actual label if different
ground_channel = find(strcmp(chanlabels(D_original), 'F9')); % Adjust 'GND' to actual label if different

% Rename reference and ground channels
if ~isempty(ref_channel)
    D_original = chanlabels(D_original, ref_channel, 'REF');
    disp(['Renamed reference channel at index ' num2str(ref_channel) ' to REF']);
end

if ~isempty(ground_channel)
    D_original = chanlabels(D_original, ground_channel, 'GND');
    disp(['Renamed ground channel at index ' num2str(ground_channel) ' to GND']);
end




%
%% STEP 3: Remove all channels except the first 3, accel, ref, and ground
%

% Identify the channels to keep
channels_to_keep = [1, 2, 3];

% Combine all channels to keep
channels_to_keep = unique([channels_to_keep, accel_channels, ref_channel, ground_channel]);

% Ensure indices are within the valid range
channels_to_keep = channels_to_keep(channels_to_keep <= D_original.nchannels);

% Create a new meeg object with only selected channels
S_crop = struct();
S_crop.D = D_original;
S_crop.channels = D_original.chanlabels(channels_to_keep);
D_cropped = spm_eeg_crop(S_crop);

% Update channels_to_keep to match the new dataset
channels_to_keep = 1:D_cropped.nchannels; % After cropping, indices reset to 1:N
disp(['Number of remaining channels: ' num2str(D_cropped.nchannels)]);
disp('Remaining channel labels:');
disp(D_cropped.chanlabels);




%
%% STEP 4: Define the montage
%

% Define montage matrix
montage_matrix = eye(D_cropped.nchannels);

% Explanation:
% if D.nchannels = 4
% montage_matrix =
%     1     0     0     0
%     0     1     0     0
%     0     0     1     0
%     0     0     0     1

% Apply the montage to keep only the desired channels
channels_to_remove = setdiff(1:D_cropped.nchannels, channels_to_keep);

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
S_mont = struct();
S_mont.D = D_cropped;
S_mont.montage.tra = montage_matrix;
S_mont.montage.labelorg = D_cropped.chanlabels;
S_mont.montage.labelnew = D_cropped.chanlabels(channels_to_keep);
S_mont.keepothers = 0;
D_montaged = spm_eeg_montage(S_mont);

disp('Applied montage successfully.');




%
%% STEP 5: Define conditions for comparisons
%

% Define the event values we need for each condition
condition_events = { % MODIFY HERE (3/13)
    % [42, 12];  % Condition 1 / LL
    [43, 13];  % Condition 2 / LH
    % [44, 14];  % Condition 3 / HL
    % [46, 16];  % Condition 4 / HH
};

% Define the condition label
condition_labels = {"Condition 2"}; % MODIFY HERE (4/13)
% condition_labels = {"Condition 1", "Condition 2", "Condition 3", "Condition 4"};




%
%% STEP 6: Fix the triggers
%

% Time window for epochs (in milliseconds)
prestim = 0; % Time before the stimulus
poststim = 1000; % Time after the stimulus

% Initialize an empty trial structure
trials = []; % Store trial definitions
trl = []; % Store trial timing (start, end, offset)
conditionlabels = {};

% Loop through each condition to define trial events
for cond = 1:length(condition_events)
    start_event = condition_events{cond}(1);
    end_event = condition_events{cond}(2);
    label = char(condition_labels{cond});
    
    % Filter the events for the current condition’s start markers
    cond_table = event_data(strcmp(event_data.type, 'Stimulus') & (event_data.value == start_event | event_data.value == end_event), :);

    % Skip to the next condition if no relevant events are found
    if isempty(cond_table)
        continue
    end

    %
    % STEP 6.1: Define trial structure
    %

    trials = [trials,
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', char([label])),
        struct('eventtype', 'Stimulus', 'eventvalue', num2str(start_event), 'prestim', prestim, 'poststim', poststim, 'conditionlabel', char([append(label, "_TEST")]))
    ];

    %
    % STEP 6.2: Generate Labels
    %
    
    % Loop over each start event and define a trial
    for i = 1:2:height(cond_table)
        if i == height(cond_table)-1
            appendedLabel = append(label, "_TEST");
            % conditionlabels = [conditionlabels, append(label, "_TEST")];
            conditionlabels = [conditionlabels, {appendedLabel}];
        else
            conditionlabels = [conditionlabels, {label}];
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
%% STEP 7: Epoch our data into 1 second intervals (to match the engagement paper)
%

fprintf('Starting epoching process...\n');

% Validate trl matrix
if any(trl(:, 2) <= trl(:, 1))
    error('Invalid trl matrix: End times must be greater than start times.');
end

if any(trl(:) < 0) || any(trl(:, 2) > D_montaged.nsamples)
    error('Invalid trl matrix: Times must be within the data range.');
end

% Sampling frequency
fs = D_montaged.fsample;

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
%% STEP 8: Baseline Correction
%

% Apply baseline correction to each epoch
% Define a structure to hold the baseline correction parameters
% S = struct();
S_bc = [];
S_bc.D = D_montaged;
S_bc.Dbaseln = D_montaged;
S_bc.timewin = [2762.35, 2777.488];  % Set baseline window (start and end samples) % MODIFY HERE (5/13)
D_bc = spm_eeg_bc(S_bc);  % Apply baseline correction

disp('Baseline correction applied successfully.');





%
%% STEP 9: Output
%

disp('TRL:');
disp(trl);
disp('Condition Labels:');
disp(conditionlabels);
disp('Trials Structure');
disp(trials);



% Set up the SPM structure
S_epochs = struct();
S_epochs.D = D_bc;
% This was here previously but I don't know if it's doing what we want it to
trl = condition_timestamps;

fs = D_bc.fsample;
S_epochs.trl = trl;
S_epochs.trl(:, 1:2) = round(S_epochs.trl(:, 1:2) * fs / 1000);

if size(S_epochs.trl, 1) == length(conditionlabels)
    disp('True');
end

conditionlabels = cellstr(condition_labels);
size(conditionlabels)
S_epochs.conditionlabels = conditionlabels;

% Background info
% Previously we had a warning saying: "There was no baseline specified. The
% data is not baseline-correct." Since I am using spm_eeg_bc above, we are
% successfully supplying the baseline correction, but spm_eeg_epochs has
% its own optional baseline. I turned off the optional baseline here.
S_epochs.bc = 0;

modified_D = spm_eeg_epochs(S_epochs);
modified_D.save();

modified_D_2 = modified_D; % MODIFY HERE (6/13)
trl_2 = trl; % MODIFY HERE (7/13)



% Set up the SPM structure for epoching
epoch_S = struct();

epoch_S.D = D_bc;
epoch_S.trl = epoch_trl;
epoch_S.conditionlabels = conditionlabels;

epoch_D = spm_eeg_epochs(epoch_S);

epoch_D.save();
epoch_D_2 = epoch_D; % MODIFY HERE (8/13)
epoch_trl_2 = epoch_trl; % MODIFY HERE (9/13)



%
%% STEP 10: Save Files in New Folder
%

% Define the parent directory for your project
preparedDataDir = 'C:/Users/Jalynn/Documents/GitHub/EEG_SPM_Analysis/Prepared Data/P11';

save(fullfile(preparedDataDir, 'modified_D_2.mat'), 'modified_D_2'); % MODIFY HERE (10/13)

writematrix(trl_2, fullfile(preparedDataDir, 'trl_2.csv')); % MODIFY HERE (11/13)

save(fullfile(preparedDataDir, 'epoch_D_2.mat'), 'epoch_D_2'); % MODIFY HERE (12/13)

writematrix(epoch_trl_2, fullfile(preparedDataDir, 'epoch_trl_2.csv')); % MODIFY HERE (13/13)

disp('Files saved successfully in the Prepared Data folder.');
