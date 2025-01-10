%
%% PREPARE.M
%

format shortG




%
%% MODIFY HERE
%

base_path = 'C:/Users/Jalynn/Desktop/BDF_EEG_Files/P14';
spm_file = fullfile(base_path, 'spmeeg_P14BDF.mat');
event_file = fullfile(base_path, 'P14events_data.csv');
condition_number = 2; % Either 1 or 2
baseline_window = [1606.594, 1621.704]; % Will automate soon
preparedDataDir = 'C:/Users/Jalynn/Documents/GitHub/EEG_SPM_Analysis/Prepared Data/P14';




%
%% STEP 0: Convert from .XDF
%

% Recall that the data initially came in an .xdf format
% This format was not automatically accepted by SPM12
% It was first converted to .mat and .dat with EEGLab
% Since we attempted so many different file conversions...
% The load function below uses spmeeg_PXXBDF.mat
% So ignore the BDF but the .mat is important




%
%% STEP 1: Load the EEG file & event data CSV
%

% Load the SPM EEG file
D_original = spm_eeg_load(spm_file);

% Load the events from CSVs
event_data = readtable(event_file);





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
if condition_number == 1
    condition_events = {[42, 12]};
    disp('Using events: [42, 12]'); % Print for condition 1
else
    condition_events = {[43, 13]};
    disp('Using events: [43, 13]'); % Print for other conditions
end

% Define the condition label
% condition_labels = {sprintf("Condition %d", condition_number)};




%
%% STEP 6: Fix the triggers
%


%
% STEP 6.1: Create table of only start and end markers
%


start_event = condition_events{1}(1); % 43 or 42
end_event = condition_events{1}(2); % 13 or 12
% label = char(condition_labels{1}); % Condition 1/2

% Filter the events for the current condition’s - only 4X and 1X values
cond_table = event_data(strcmp(event_data.type, 'Stimulus') & (event_data.value == start_event | event_data.value == end_event), :);

% Print error if no relevant events are found
if isempty(cond_table)
    error('Condition table is empty. Please check the input data.');
end


%
% STEP 6.2: Generate TRL & TRL Labels with 7 Trials
%


trl = []; % Store 7 trials (start time, end time, offset)
trl_labels = {}; % Store labels for 7 trials (Condition X (x6) & Condition X_Test)

% Loop over each start event and define a trial
for i = 1:2:height(cond_table)
    if i == height(cond_table)-1
        trl_labels = [trl_labels, {sprintf("Condition %d_TEST", condition_number)}];
    else
        trl_labels = [trl_labels, {sprintf("Condition %d", condition_number)}];
    end

    % Store trial start, end, and offset
    trl = [trl; cond_table{i, "time"} cond_table{i+1, "time"} 0];
end


%
% STEP 6.3: Extract Start & End of Train/Test Timestamps
%


%
% STEP 6.3.1: Training Start & End Time
%


% Initialize a structure to store the results
training_timestamps = [];

% Find the first occurrence of the start event
first_start_idx = find(event_data.value == start_event, 1, 'first');

% Find the second to last occurence of the end event
all_end_idx = find(event_data.value == end_event);
if numel(all_end_idx) >= 2
    second_last_end_idx = all_end_idx(end - 1);
end

% Save timestamps if both indices are valid
if ~isempty(first_start_idx) && ~isempty(second_last_end_idx)
    first_start_time = event_data.time(first_start_idx);
    second_last_end_time = event_data.time(second_last_end_idx);
    
    % Append to condition_timestamps
    training_timestamps = [training_timestamps; first_start_time, second_last_end_time, 0];
else
    fprintf('Event indices not found');
end


%
% STEP 6.3.2: Testing Start & End Time
%


% Initialize a structure to store the results
testing_timestamps = [];

% Find the last occurent of the start event
last_start_idx = find(event_data.value == start_event, 1, 'last');

% Find the last occurence of the end event
last_end_idx = find(event_data.value == end_event, 1, 'last');

% Save timestamps if both indices are valid
if ~isempty(second_last_start_idx) && ~isempty(last_end_idx)
    last_start_time = event_data.time(last_start_idx);
    last_end_time = event_data.time(last_end_idx);
    
    % Append to condition_timestamps
    testing_timestamps = [testing_timestamps; last_start_time, last_end_time, 0];
else
    fprintf('Event indices not found');
end


%
% STEP 6.3.3: Total Condition Start & End Time
%


% Initialize a structure to store the results
condition_timestamps = [];

% Save timestamps if both indices are valid
if ~isempty(first_start_idx) && ~isempty(last_end_idx)
    first_start_time = event_data.time(first_start_idx);
    last_end_time = event_data.time(last_end_idx);
    
    % Append to condition_timestamps
    condition_timestamps = [condition_timestamps; first_start_time, last_end_time, 0];
else
    fprintf('Event indices not found');
end




%
%% STEP 7: Epoch our data into 1 second intervals (to match the engagement paper)
%

fprintf('Starting epoching process...\n');

% Validate 7 trial trl matrix
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

% new trl for 1 second epoching - All 7 trials
epoch_trl = [];

% Generate 1 second epoch intervals
for i = 1:size(trl, 1)
    start_time = trl(i, 1);
    end_time = trl(i, 2);
    offset = trl(i, 3);

    % Create 1-second intervals between start_time and end_time
    for t = start_time:1:end_time - 1  % Increment by 1 second, ensuring we don't exceed the end_time
        epoch_trl = [epoch_trl; t, t + 1, offset];
    end
end

% Display a message if epoch_trl is still empty after processing
if isempty(epoch_trl)
    warning('No valid 1-second intervals were created. Verify trl data.');
else
    fprintf('Epoching completed successfully. Total epochs: %d\n', size(epoch_trl, 1));
end




%
%% STEP 8: Baseline Correction
%

% Apply baseline correction to each epoch
% Define a structure to hold the baseline correction parameters
S_bc = [];
S_bc.D = D_montaged;
S_bc.Dbaseln = D_montaged;
S_bc.timewin = baseline_window; % Set baseline window (start and end samples)
D_bc = spm_eeg_bc(S_bc);  % Apply baseline correction

disp('Baseline correction applied successfully.');





%
%% STEP 9: Output
%


%
% STEP 9.1: Command Window Printing
%


% CSV File was broken into one condition table

disp('Condition Table:');
disp(cond_table);

% There is only one complete trl for all 7 trials

disp('7 Trial TRL:');
disp(trl);

% There is one corresponding trl_labels for the 7 trial labels

disp('7 Trial Labels:');
disp(trl_labels);

% There is only one complete epoch_trl for 1s intervals for all 7 trials

disp('Epoch TRL:');
disp(epoch_trl);

% All 7 trials in one start/end line called condition_timestamps

disp('Condition Timestamps:');
disp(condition_timestamps);

% First 6 trails in one start/end line called training_timestamps

disp('Training Timestamps:');
disp(training_timestamps);

% Last trial in one start/end line called testing_timestamps

disp('Testing Timestamps:');
disp(testing_timestamps);

% TRL / TRL_Labels & Epoch_TRL still needs train/test splitting
% D_bc still needs train/test splitting


%
% STEP 9.2: continuous_D
%

% Not ideal but working
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

% added since the updated script removed it above
condition_labels = {sprintf("Condition %d", condition_number)};

conditionlabels = cellstr(condition_labels);
size(conditionlabels)
S_epochs.conditionlabels = conditionlabels;

S_epochs.bc = 0;

modified_D = spm_eeg_epochs(S_epochs);
modified_D.save();

modified_D_2 = modified_D; % MODIFY HERE (6/14)
trl_2 = trl; % MODIFY HERE (7/14)


% cont_S = struct();
% cont_S.D = D_bc;
% 
% % Recall TRL is 7 lines of start and end times
% % We can't use this because they are not equal times
% % Complete (not train/test) data uses the complete condition timestamps
% cont_S.trl = condition_timestamps;
% eval(sprintf('trl_%d = condition_timestamps;', condition_number));
% 
% % Not sure what this was used for yet
% % fs = D_bc.fsample;
% % cont_S.trl(:, 1:2) = round(cont_S.trl(:, 1:2) * fs / 1000);
% % if size(cont_S.trl, 1) == length(trl_labels)
% %     disp('True');
% % end
% 
% % Removed condition_labels and conditionlabels
% % trl_labels = cellstr(condition_labels);
% % size(trl_labels)
% 
% % One line trl - one line condition label
% condition_label = trl_labels{1};
% cont_S.conditionlabels = condition_label; 
% 
% % Background info
% % Previously we had a warning saying: "There was no baseline specified. The
% % data is not baseline-correct." Since I am using spm_eeg_bc above, we are
% % successfully supplying the baseline correction, but spm_eeg_epochs has
% % its own optional baseline. I turned off the optional baseline here.
% cont_S.bc = 0;
% 
% continuous_D = spm_eeg_epochs(cont_S);
% continuous_D.save();
% 
% eval(sprintf('continuous_D_%d = continuous_D;', condition_number))




%
% STEP 9.2.1: train_continuous_D
%


%
% STEP 9.2.2: test_continuous_D
%


%
% STEP 9.3: epoch_D
%

% not ideal but working

% Set up the SPM structure for epoching
epoch_S = struct();

epoch_S.D = D_bc;
epoch_S.trl = epoch_trl;
epoch_S.conditionlabels = conditionlabels;

epoch_D = spm_eeg_epochs(epoch_S);

epoch_D.save();
epoch_D_2 = epoch_D; % MODIFY HERE (8/14)
epoch_trl_2 = epoch_trl; % MODIFY HERE (9/14)

% epoch_S = struct();
% 
% epoch_S.D = D_bc;
% epoch_S.trl = epoch_trl;
% epoch_S.conditionlabels = trl_labels;
% 
% epoch_D = spm_eeg_epochs(epoch_S);
% 
% epoch_D.save();
% eval(sprintf('epoch_D_%d = epoch_D;', condition_number));
% eval(sprintf('epoch_trl_%d = epoch_trl;', condition_number));


%
% STEP 9.3.1 train_epoch_D
%


%
% STEP 9.3.2: test_epoch_D
%




%
%% STEP 10: Save Files in New Folder
%

% Define the parent directory for your project
save(fullfile(preparedDataDir, sprintf('continuous_D_%d.mat', condition_number)), sprintf('continuous_D_%d', condition_number));

writematrix(eval(sprintf('trl_%d', condition_number)), fullfile(preparedDataDir, sprintf('trl_%d.csv', condition_number)));

save(fullfile(preparedDataDir, sprintf('epoch_D_%d.mat', condition_number)), sprintf('epoch_D_%d', condition_number));

writematrix(eval(sprintf('epoch_trl_%d', condition_number)), fullfile(preparedDataDir, sprintf('epoch_trl_%d.csv', condition_number)));

disp('Files saved successfully in the Prepared Data folder.');
