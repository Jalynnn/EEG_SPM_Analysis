D = spm_eeg_load('spmeeg_01.mat');

% Retrieve the structure of the existing events
current_events = D.events;

new_events = struct('type', {}, 'value', {}, 'time', {}, 'duration', {}, 'offset', {});

% Example: Adding start and end markers for each condition
% Replace <time_of_event> with actual time in seconds

% Condition 1
new_events(end+1) = struct('type', 'LSL', 'value', 46, 'time', 680863, 'duration', 0, 'offset', 0);
new_events(end+1) = struct('type', 'LSL', 'value', 16, 'time', 954886, 'duration', 0, 'offset', 0);

% Condition 2
new_events(end+1) = struct('type', 'LSL', 'value', 42, 'time', 1086757, 'duration', 0, 'offset', 0);
new_events(end+1) = struct('type', 'LSL', 'value', 12, 'time', 1260746, 'duration', 0, 'offset', 0);

% Condition 3
new_events(end+1) = struct('type', 'LSL', 'value', 43, 'time', 1367694, 'duration', 0, 'offset', 0);
new_events(end+1) = struct('type', 'LSL', 'value', 13, 'time', 1543537, 'duration', 0, 'offset', 0);

% Condition 4
new_events(end+1) = struct('type', 'LSL', 'value', 44, 'time', 1655958, 'duration', 0, 'offset', 0);
new_events(end+1) = struct('type', 'LSL', 'value', 14, 'time', 1850801, 'duration', 0, 'offset', 0);

% Not allowed to modiify d.events directly
% D.events = [D.events, new_events];

% spm_eeg_add_events doesn't exist here
% Loop through each event and add it to D using spm_eeg_add_events
% for i = 1:numel(new_events)
%     D = spm_eeg_add_events(D, new_events(i));
% end

% Again, not allowed to modify D.events
% Ensure the fields are ordered the same way
% new_events = orderfields(new_events, current_events);
% Concatenate current events with new events
% D.events = [current_events, new_events];  % Use semicolon to concatenate structs vertically

D.events = D.events;

save(D);
