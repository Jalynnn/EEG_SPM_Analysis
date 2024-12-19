% Step 1: Load the SPM EEG file with epochs already defined
%D = spm_eeg_load('/spmeeg_P07BDF.mat'); % Load the newly epoch-processed file

%% Step 2: Downsample to 250 Hz (if not already downsampled)
D = modified_D;
disp(D)

if D.fsample > 250
    S = struct('D', D, 'fsample_new', 250);
    D = spm_eeg_downsample(S);
end

%% Step 3: Apply Bandpass Filter (0.1 Hz to 30 Hz)
S = struct('D', D, 'band', 'bandpass', 'freq', [0.1 30]);
D = spm_eeg_filter(S);

%% Step 4: Artifact Rejection on the Defined Epochs
% Define methods for artifact rejection
S = struct();
S.D = D;
S.methods(1).channels = 'EEG';           % Apply to EEG channels only
S.methods(1).fun = 'flat';               % Flatline detection
S.methods(1).settings.threshold = 0.001; % Flatline threshold in µV
S.methods(1).settings.seqlength = 10;    % Minimum sequence length in samples for a flatline

S.methods(2).channels = 'EEG';           % Apply to EEG channels only
S.methods(2).fun = 'peak2peak';          % Peak-to-peak thresholding
S.methods(2).settings.threshold = 200;   % Peak-to-peak threshold in µV

% Apply artifact rejection using the defined methods
D = spm_eeg_artefact(S);

%% Step 5: Calculate Power in Alpha and Theta Bands for Each Epoch
% Define the conditions to calculate ratios for
condition_labels = {'Condition 1'};

% Frequency bands for alpha and theta
alpha_band = [8 13];
theta_band = [4 7];

% Sampling frequency
fs = D.fsample;

% Initialize an array to store alpha-theta ratios for each condition
alpha_theta_ratios = zeros(1, length(condition_labels));

% Loop over each condition to calculate the alpha-theta ratio
for cond = 1:length(condition_labels)
    % Find epochs that match the current condition label
    epoch_indices = find(strcmp(D.conditions, condition_labels{cond}));
    
    % Check if there are any matching epochs for the condition
    if isempty(epoch_indices)
        fprintf('No epochs found for %s. Setting alpha-theta ratio to NaN.\n', condition_labels{cond});
        alpha_theta_ratios(cond) = NaN;
        continue;
    end
    
    % Preallocate arrays to store power values for alpha and theta bands
    alpha_power = zeros(1, length(epoch_indices));
    theta_power = zeros(1, length(epoch_indices));

    % Calculate power for each epoch in the current condition
    for e = 1:length(epoch_indices)
        % Get the data for the current epoch
        epoch_data = D(:, :, epoch_indices(e)); % Extract data for the epoch

        % Check for empty or NaN data in the epoch
        if all(isnan(epoch_data(:)))
            fprintf('Epoch %d for %s contains NaNs or is empty.\n', epoch_indices(e), condition_labels{cond});
            alpha_power(e) = NaN;
            theta_power(e) = NaN;
            continue;
        end

        % Calculate power spectral density for the epoch
        [pxx, f] = pwelch(epoch_data', [], [], [], fs);

        % Find indices for alpha and theta bands
        alpha_idx = (f >= alpha_band(1)) & (f <= alpha_band(2));
        theta_idx = (f >= theta_band(1)) & (f <= theta_band(2));

        % Calculate mean power in each band for the epoch
        alpha_power(e) = mean(pxx(alpha_idx, :), 'all'); % Mean alpha power across channels
        theta_power(e) = mean(pxx(theta_idx, :), 'all'); % Mean theta power across channels
    end

    % Remove NaN values before calculating the final ratio
    alpha_power = alpha_power(~isnan(alpha_power));
    theta_power = theta_power(~isnan(theta_power));

    % Calculate the alpha-theta ratio for the current condition, if valid data exists
    if ~isempty(alpha_power) && ~isempty(theta_power) && mean(theta_power) > 0
        alpha_theta_ratios(cond) = mean(alpha_power) / mean(theta_power);
    else
        alpha_theta_ratios(cond) = NaN; % Set NaN if no valid data exists
    end
end

% Display the alpha-theta ratio for each condition
for cond = 1:length(condition_labels)
    fprintf('%s Alpha-Theta Ratio: %.3f\n', condition_labels{cond}, alpha_theta_ratios(cond));
end