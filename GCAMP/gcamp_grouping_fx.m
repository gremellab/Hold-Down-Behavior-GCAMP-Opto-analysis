
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
GCAMP_Save_Dir = 'C:\Users\dcschrei\Desktop\laptop transfer\0s gone regression data and excel\gcampdata\m2 reshuffle';
cd(GCAMP_Save_Dir)
indiv_gcamps = dir('*.mat');
num_gcamps = numel(indiv_gcamps);

for i = 1:num_gcamps
    load([indiv_gcamps(i).folder '\' indiv_gcamps(i).name])
    Grouped_GCAMP.Mice{i} = GCAMP;
end

% save([GCAMP_Save_Dir 'GCAMP_Grouped-' Grouped_GCAMP.Mice{1}.training_day], 'GCAMP');
save(['GCAMP_Grouped-' Grouped_GCAMP.Mice{1}.training_day], 'Grouped_GCAMP','-v7.3');

base_time_start = -15;
base_time_end = -5;
time_end = 5;
GCAMP_plot_with_baseline_Z_scores(Grouped_GCAMP, base_time_start, base_time_end,time_end)
% GCAMP_plot_with_baseline_Z_scores_prob(Grouped_GCAMP, base_time_start, base_time_end,time_end)

%  GCAMP_grand_regression(Grouped_GCAMP, base_time_start, base_time_end,time_end)

