%
%% Step 1: Load the SPM EEG file with epochs already defined
%

% This was completed in prepare.m
% D = spm_eeg_load('/spmeeg_P07BDF.mat'); % Load the newly epoch-processed file

% MODIFY HERE (1/4)
D = epoch_D_2; % continue without epochs for now
% D = test_continuous_D_2;

disp(D)




%
%% Step 2: Downsample to 250 Hz (if not already downsampled)
%

if D.fsample > 250
    S = struct('D', D, 'fsample_new', 250);
    D = spm_eeg_downsample(S);
end




%
%% Step 3: Apply Bandpass Filter (0.1 Hz to 30 Hz)
%

S = struct('D', D, 'band', 'bandpass', 'freq', [0.1 30]);
D = spm_eeg_filter(S);




%
%% Step 4: Artifact Rejection on the Defined Epochs
%

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




%
%% Step 5: Calculate Power in the various bands for ratios
% the alpha/theta uses specific channels (Pz-Alpha and Fz-Theta)
% the beta/(a+t) is calc'd at every channel and summed
%

% Define the conditions to calculate ratios for
condition_labels = {'Condition 2'}; % MODIFY HERE (2/4)

% Frequency bands for alpha, beta and theta
beta_band = [14 20];
alpha_band = [8 13];
theta_band = [4 7];

% continuous data 

% Frequency bands for alpha, beta and theta
beta_band = [14 20];
alpha_band = [8 13];
theta_band = [4 7];

% Sampling frequency
fs = D.fsample;

% Get the data for all channels (continuous)
data = D(:,:);

% Initialize variables for storing power values
cont_pz_alpha_power = [];
cont_fz_theta_power = [];
cont_beta_power = [];
cont_alpha_power = [];
cont_theta_power = [];

% Calculate power spectral density (PSD) for the continuous data
[pxx, f] = pwelch(data', [], [], [], fs);

% Find indices corresponding to the defined frequency bands
cont_alpha_idx = (f >= alpha_band(1)) & (f <= alpha_band(2));
cont_theta_idx = (f >= theta_band(1)) & (f <= theta_band(2));
cont_beta_idx = (f >= beta_band(1)) & (f <= beta_band(2));

% Compute power for specific channels and all channels
if any(strcmp(D.chanlabels, 'Pz')) && any(strcmp(D.chanlabels, 'Fz'))
    % Alpha power for Pz channel
    cont_pz_alpha_power = mean(pxx(cont_alpha_idx, strcmp(D.chanlabels, 'Pz')));
    % Theta power for Fz channel
    cont_fz_theta_power = mean(pxx(cont_theta_idx, strcmp(D.chanlabels, 'Fz')));
else
    error('Pz or Fz channel not found in the data.');
end

% Compute mean power for beta, alpha, and theta across all channels
cont_beta_power = mean(pxx(cont_beta_idx, :), 'all'); 
cont_alpha_power = mean(pxx(cont_alpha_idx, :), 'all'); 
cont_theta_power = mean(pxx(cont_theta_idx, :), 'all'); 

% Calculate the Pz Alpha / Fz Theta ratio
if ~isempty(cont_pz_alpha_power) && ~isempty(cont_fz_theta_power) && cont_fz_theta_power > 0
    cont_pz_alpha_fz_theta_ratio = cont_pz_alpha_power / cont_fz_theta_power;
else
    cont_pz_alpha_fz_theta_ratio = NaN;
end

% Calculate Beta / (Alpha + Theta) ratio averaged across all channels
if ~isempty(cont_beta_power) && ~isempty(cont_alpha_power) && ~isempty(cont_theta_power) && (cont_alpha_power + cont_theta_power > 0)
    cont_beta_over_alpha_theta_ratio = cont_beta_power / (cont_alpha_power + cont_theta_power);
else
    cont_beta_over_alpha_theta_ratio = NaN;
end

% Display results
fprintf('Pz Alpha / Fz Theta Ratio: %.3f\n', cont_pz_alpha_fz_theta_ratio);
fprintf('Beta / (Alpha + Theta) Ratio (All Channels): %.3f\n', cont_beta_over_alpha_theta_ratio);

% Save results to a file
output_file = fullfile(preparedDataDir, 'continuous_alpha_theta_beta_ratios_output.txt');
fileID = fopen(output_file, 'w');
fprintf(fileID, 'Pz Alpha / Fz Theta Ratio: %.3f\n', cont_pz_alpha_fz_theta_ratio);
fprintf(fileID, 'Beta / (Alpha + Theta) Ratio (All Channels): %.3f\n', cont_beta_over_alpha_theta_ratio);
fclose(fileID);


% epoched data 

% Sampling frequency
fs = D.fsample;

% Initialize variables for storing results
pz_alpha_fz_theta_ratios = zeros(1, length(condition_labels));
beta_over_alpha_theta_ratios = zeros(1, length(condition_labels));

% Loop over each condition
for cond = 1:length(condition_labels)
    % Find epochs that match the current condition label
    epoch_indices = find(strcmp(D.conditions, condition_labels{cond}));

    if isempty(epoch_indices)
        fprintf('No epochs found for %s. Setting ratios to NaN.\n', condition_labels{cond});
        pz_alpha_fz_theta_ratios(cond) = NaN;
        beta_over_alpha_theta_ratios(cond) = NaN;
        continue;
    end

    % Initialize storage for Pz alpha, Fz theta, beta, alpha, and theta power
    ep_pz_alpha_power = [];
    ep_fz_theta_power = [];
    ep_beta_power = [];
    ep_alpha_power = [];
    ep_theta_power = [];

    % Process each epoch for the current condition
    for e = 1:length(epoch_indices)
        epoch_data = D(:, :, epoch_indices(e)); % Get epoch data

        if all(isnan(epoch_data(:)))
            continue;
        end

        % Calculate power spectral density (PSD)
        [pxx, f] = pwelch(epoch_data', [], [], [], fs);

        % Extract indices corresponding to the defined frequency bands
        ep_alpha_idx = (f >= alpha_band(1)) & (f <= alpha_band(2));
        ep_theta_idx = (f >= theta_band(1)) & (f <= theta_band(2));
        ep_beta_idx = (f >= beta_band(1)) & (f <= beta_band(2));

        % Calculate power for specific channels and all channels
        if any(strcmp(D.chanlabels, 'Pz')) && any(strcmp(D.chanlabels, 'Fz'))
            % Alpha power for Pz channel
            ep_pz_alpha_power = [ep_pz_alpha_power, mean(pxx(ep_alpha_idx, strcmp(D.chanlabels, 'Pz')))];
            % Theta power for Fz channel
            ep_fz_theta_power = [ep_fz_theta_power, mean(pxx(ep_theta_idx, strcmp(D.chanlabels, 'Fz')))];
        end

        % Compute power for beta, alpha, and theta across all channels
        ep_beta_power = [ep_beta_power, mean(pxx(ep_beta_idx, :), 'all')];
        ep_alpha_power = [ep_alpha_power, mean(pxx(ep_alpha_idx, :), 'all')];
        ep_theta_power = [ep_theta_power, mean(pxx(ep_theta_idx, :), 'all')];
    end

    % Calculate Pz Alpha / Fz Theta ratio
    if ~isempty(ep_pz_alpha_power) && ~isempty(ep_fz_theta_power) && mean(ep_fz_theta_power) > 0
        pz_alpha_fz_theta_ratios(cond) = mean(ep_pz_alpha_power) / mean(ep_fz_theta_power);
    else
        pz_alpha_fz_theta_ratios(cond) = NaN;
    end

    % Calculate Beta / (Alpha + Theta) ratio averaged across all channels
    if ~isempty(ep_beta_power) && ~isempty(ep_alpha_power) && ~isempty(ep_theta_power) && (mean(ep_alpha_power) + mean(ep_theta_power) > 0)
        beta_over_alpha_theta_ratios(cond) = mean(ep_beta_power) / (mean(ep_alpha_power) + mean(ep_theta_power));
    else
        beta_over_alpha_theta_ratios(cond) = NaN;
    end
end

% Display results
for cond = 1:length(condition_labels)
    fprintf('%s Pz Alpha / Fz Theta Ratio: %.3f\n', condition_labels{cond}, pz_alpha_fz_theta_ratios(cond));
    fprintf('%s Beta / (Alpha + Theta) Ratio (All Channels): %.3f\n', condition_labels{cond}, beta_over_alpha_theta_ratios(cond));
end

% Save results to a file
output_file = fullfile(preparedDataDir, 'epoched_alpha_theta_beta_ratios_output.txt');
fileID = fopen(output_file, 'w');
for cond = 1:length(condition_labels)
    fprintf(fileID, '%s Pz Alpha / Fz Theta Ratio: %.3f\n', condition_labels{cond}, pz_alpha_fz_theta_ratios(cond));
    fprintf(fileID, '%s Beta / (Alpha + Theta) Ratio (All Channels): %.3f\n', condition_labels{cond}, beta_over_alpha_theta_ratios(cond));
end
fclose(fileID);






% Save to a text file

% MODIFY HERE (3/4)
% Specify the directory where the file should be saved
preparedDataDir = 'C:/Users/Jalynn/Documents/GitHub/EEG_SPM_Analysis/Prepared Data/P14';

% MODIFY HERE (4/4)
% Combine the directory path with the desired file name
file_path = fullfile(preparedDataDir, 'MD2_alpha_theta_ratios_output.txt');

% Open the file for writing
fileID = fopen(file_path, 'w');

% % Display the alpha-theta ratio for each condition and write to the file
% for cond = 1:length(condition_labels)
%     fprintf(fileID, '%s Alpha-Theta Ratio: %.3f\n', condition_labels{cond}, alpha_theta_ratios(cond));
% end

% Close the file
fclose(fileID);


