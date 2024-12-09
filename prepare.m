format shortG

% Load the SPM EEG file
D = spm_eeg_load('P10/spmeeg_P10BDF.mat');

% Load the events from CSV
event_data = readtable('P10/P10events_data.csv'); % Replace with the actual path

% Define the event values we need for each condition
condition_events = {
    [42, 12];  % Condition 1
    [43, 13];  % Condition 2
    [44, 14];  % Condition 3
    [46, 16];  % Condition 4
};
condition_labels = {"Condition 1", "Condition 2", "Condition 3", "Condition 4"};

% Sampling frequency
fs = D.fsample;

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


trl
conditionlabels
trials

trl = condition_timestamps
remove all conditions except the 42/12

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
