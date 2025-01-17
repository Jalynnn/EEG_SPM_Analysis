% I recreated this script to only take in the two conditions and modified
% the PSD heavily since we have two D variables instead of one.

%
%% Step 1: Load the epoched EEG data file with labeled conditions
%

condition_labels = {'Condition 1', 'Condition 2'};

% Sampling frequency
% fs = D.fsample;
fs = 250;

% Define frequency bands of interest
freq_bands = struct('theta', [4 7], 'alpha', [8 13], 'beta', [13 30]);
freq_band_names = fieldnames(freq_bands);




%
%% Step 2: Power Spectral Density (PSD) Plot
%

figure;
subplot(1, 2, 1); % Create a 1x2 subplot layout
hold on;

max_frequency = 50; % Limit frequency range to 50 Hz for better visualization

% for cond = 1:length(condition_labels)
%     % Find epochs that match the current condition label
%     % epoch_indices = find(strcmp(D.conditions, condition_labels{cond}));
%     % current_data = condition_data{cond};
% 
%     % % Check if there are any matching epochs for the condition
%     % if isempty(epoch_indices)
%     % % if isempty(current_data)
%     %     fprintf('No epochs found for %s. Skipping PSD calculation.\n', condition_labels{cond});
%     %     continue;
%     % end
% 
%     % Extract data for the current condition based on the identified epochs
%     % condition_data = D(:, :, epoch_indices);
%     condition_data
% 
%     % Calculate PSD using MATLAB's pwelch
%     [psd, f] = pwelch(mean(condition_data, 3)', [], [], [], fs);
%     % [psd, f] = pwelch(mean(current_data, 3)', [], [], [], fs);
% 
%     % Limit the frequency range for plotting
%     freq_range = f <= max_frequency;
% 
%     % Plot the PSD
%     plot(f(freq_range), mean(psd(freq_range, :), 2), 'DisplayName', condition_labels{cond});
% end

condition_data_1 = modified_D_1(:, :, :);
[psd, f] = pwelch(mean(condition_data_1, 3)', [], [], [], fs);
% Limit the frequency range for plotting
freq_range = f <= max_frequency;
% Plot the PSD
plot(f(freq_range), mean(psd(freq_range, :), 2), 'DisplayName', condition_labels{1});

condition_data_2 = modified_D_2(:, :, :);
[psd, f] = pwelch(mean(condition_data_2, 3)', [], [], [], fs);
% Limit the frequency range for plotting
freq_range = f <= max_frequency;
% Plot the PSD
plot(f(freq_range), mean(psd(freq_range, :), 2), 'DisplayName', condition_labels{2});

hold off;
legend(condition_labels); % Use the correct condition labels for the legend
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectral Density per Condition');




%
%% Step 3: Band-Specific Power Bar Plot
% 

subplot(1, 2, 2); % Plot the bar chart in the second subplot
band_powers = zeros(length(condition_labels), length(freq_band_names));

% for cond = 1:length(condition_labels)
%     epoch_indices = find(strcmp(D.conditions, condition_labels{cond}));
% 
%     % Skip if no epochs found for the condition
%     if isempty(epoch_indices)
%         fprintf('No epochs found for %s. Skipping band-specific power calculation.\n', condition_labels{cond});
%         band_powers(cond, :) = NaN;
%         continue;
%     end

for b = 1:length(freq_band_names)
    band = freq_band_names{b};
    freq_range = freq_bands.(band);
    
    % % Calculate mean band power for each condition
    % condition_data = mean(D(:, :, epoch_indices), 3);
    % filtered_data = bandpass(condition_data', freq_range, fs);
    % band_powers(cond, b) = mean(mean(filtered_data.^2, 1), 'all'); % Mean power across channels

    % Calculate mean band power for each condition
    condition_data = mean(modified_D_1(:, :, :), 3);
    filtered_data = bandpass(condition_data', freq_range, fs);
    band_powers(cond, b) = mean(mean(filtered_data.^2, 1), 'all'); % Mean power across channels

    % Calculate mean band power for each condition
    condition_data = mean(modified_D_2(:, :, :), 3);
    filtered_data = bandpass(condition_data', freq_range, fs);
    band_powers(cond, b) = mean(mean(filtered_data.^2, 1), 'all'); % Mean power across channels
end
% end

% Create the bar plot for band-specific power
bar(band_powers);
set(gca, 'xticklabel', condition_labels); % Use correct condition labels
legend(freq_band_names, 'Location', 'NorthEast');
xlabel('Condition');
ylabel('Power');
title('Band-Specific Power per Condition');
