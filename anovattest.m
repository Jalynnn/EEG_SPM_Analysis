% Sample data: Alpha-Theta Ratios for all participants across conditions
% Replace these arrays with your actual participant data for both conditions.

% P11, P12, P13, & P14
condition_1 = [0.080, 0.609, 0.531, 0.348];  % Condition 1 Alpha-Theta Ratios (all participants)
condition_2 = [0.160, 0.454, 0.492, 0.382];  % Condition 2 Alpha-Theta Ratios (all participants)

% Step 1: Perform a One-Way Repeated Measures ANOVA
% Create a table with the alpha-theta ratios for both conditions
participant_ids = (1:length(condition_1))';  % Create participant ID list
alpha_theta_data = table(participant_ids, condition_1', condition_2', 'VariableNames', {'ParticipantID', 'Condition1', 'Condition2'});

% Fit a repeated measures model
rm = fitrm(alpha_theta_data, 'Condition1-Condition2 ~ 1', 'WithinDesign', [1 2]);

% Perform repeated measures ANOVA
anova_result = ranova(rm);
disp('Repeated Measures ANOVA Results:');
disp(anova_result);

% Step 2: Perform Paired t-test across all participants
% If the data is normally distributed, use a paired t-test
[h, p_value, ci, stats] = ttest(condition_1, condition_2);
disp(['Paired t-test result: t = ', num2str(stats.tstat), ', p = ', num2str(p_value)]);

% Step 3: Visualize the Alpha-Theta Ratios across conditions for all participants
figure;
hold on;
x = 1:length(condition_1);  % Participant indices
bar_width = 0.35;           % Bar width

% Plot bars for Condition 1 and Condition 2
bar(x - bar_width/2, condition_1, bar_width, 'FaceColor', 'b', 'DisplayName', 'Condition 1');
bar(x + bar_width/2, condition_2, bar_width, 'FaceColor', 'g', 'DisplayName', 'Condition 2');

% Labeling the plot
xlabel('Participants');
ylabel('Alpha-Theta Ratio');
title('Comparison of Alpha-Theta Ratios Across Conditions');
xticks(x);
xticklabels(cellstr(num2str((1:length(condition_1))')));  % Display participant numbers
legend('show');

% Formatting the plot
hold off;
grid on;
