function [] = Grouped_Vid_Data_Function_r(Grouped_Vid_Data)
%last updated 5/12/21 Drew Schreienr
%Take the grouped_Vid_data structure as input
%Generate a table with all the mouse/session data
%use that table to build LMEs to predict duration given prior behavior

%Especially interested in looking for main effects and interactions of
%light (stim). Main effects indicate direct effect of light on press
%duration (or on subsequent press durations), while interactions with n -1
%duration or moving average indicate that light affects the use of prior
%experience.

ipi_stim_all =[];
ipi_no_stim_all =[];
ipi_stim_mean_mouse =[];
ipi_no_stim_mean_mouse = [];
regression_cell_All =[];
data_table_variables = Grouped_Vid_Data.Mice{1}.data_table_variables;

%varialbes for first press after stim or no based on time poximity
nearest_lp_latency_first_stim_mean_mouse =[];
nearest_lp_latency_first_no_stim_mean_mouse =[];
nearest_lp_latency_first_stim_all =[];
nearest_lp_latency_first_no_stim_all =[];
nearest_lp_duration_first_stim_mean_mouse =[];
nearest_lp_duration_first_no_stim_mean_mouse =[];

for i = 1:length(Grouped_Vid_Data.Mice)
    
regression_cell_All = [regression_cell_All; Grouped_Vid_Data.Mice{i}.regression_cell];
ipi_stim_all =[ipi_stim_all;  Grouped_Vid_Data.Mice{i}.ipi1_stim];
ipi_no_stim_all =[ipi_no_stim_all;  Grouped_Vid_Data.Mice{i}.ipi1_nostim];
ipi_stim_mean_mouse =[ipi_stim_mean_mouse  Grouped_Vid_Data.Mice{i}.ipi1_stim_mean];
ipi_no_stim_mean_mouse = [ipi_no_stim_mean_mouse  Grouped_Vid_Data.Mice{i}.ipi1_nostim_mean];
nearest_lp_latency_first_stim_mean_mouse =[nearest_lp_latency_first_stim_mean_mouse; Grouped_Vid_Data.Mice{i}.nearest_lp_latency_first_stim_mean];
nearest_lp_latency_first_no_stim_mean_mouse =[nearest_lp_latency_first_no_stim_mean_mouse; Grouped_Vid_Data.Mice{i}.nearest_lp_latency_first_no_stim_mean];
nearest_lp_latency_first_stim_all =[nearest_lp_latency_first_stim_all; Grouped_Vid_Data.Mice{i}.nearest_lp_latency_first_stim];
nearest_lp_latency_first_no_stim_all =[nearest_lp_latency_first_no_stim_all; Grouped_Vid_Data.Mice{i}.nearest_lp_latency_first_no_stim];
nearest_lp_duration_first_stim_mean_mouse =[nearest_lp_duration_first_stim_mean_mouse; Grouped_Vid_Data.Mice{i}.mean_stim_Durations];
nearest_lp_duration_first_no_stim_mean_mouse =[nearest_lp_duration_first_no_stim_mean_mouse; Grouped_Vid_Data.Mice{i}.mean_no_stim_Durations];

end

%import regression cell and turn to table, make sure vars are categorical
Regression_table = cell2table(regression_cell_All,'VariableNames',data_table_variables);
Regression_table.n_Stim_All = categorical(Regression_table.n_Stim_All);
Regression_table.n_minus_one_Stim_All = categorical(Regression_table.n_minus_one_Stim_All);
Regression_table.n_minus_two_Stim_All = categorical(Regression_table.n_minus_two_Stim_All);
Regression_table.n_minus_three_Stim_All = categorical(Regression_table.n_minus_three_Stim_All);
Regression_table.n_minus_four_Stim_All = categorical(Regression_table.n_minus_four_Stim_All);
Regression_table.n_minus_five_Stim_All = categorical(Regression_table.n_minus_five_Stim_All);
Regression_table.n_minus_six_Stim_All = categorical(Regression_table.n_minus_six_Stim_All);
Regression_table.mouse_indicator = categorical(Regression_table.mouse_indicator);
Regression_table.day_indicator = categorical(Regression_table.day_indicator);
Regression_table.n_minus_one_All = categorical(Regression_table.n_minus_one_All);
Regression_table.HE_n1_Indicator = categorical(Regression_table.HE_n1_Indicator);
Regression_table.stim_and_zone_indicator = categorical(Regression_table.stim_and_zone_indicator);
Regression_table.zone_indicator = categorical(Regression_table.zone_indicator);

%find any really long presses (>10s) and exclude from the model
longidx =find(Regression_table.n_duration>=10000);


 %% LME models
 
%because n_Stim is indicator whether a stim happened before press n, it
%should go into the interaction with n - 1 (question being asked: does the
%n and n - 1 relationship change if stim occurs in between). Thus all the
%interactions will appear to be "one off"

%% beh only model
modelspec_n6_only = 'n_duration ~ lp_on_times + criteria_percent_indicator +  n_minus_one_Durations_All+ n_minus_two_Durations_All + n_minus_three_Durations_All+ n_minus_four_Durations_All+ n_minus_five_Durations_All + n_minus_six_Durations_All + n_minus_seven_Durations_All + n_minus_eight_Durations_All + n_minus_nine_Durations_All + n_minus_ten_Durations_All +(1|mouse_indicator) + (1|day_indicator)';
n6_only_LME = fitlme(Regression_table, modelspec_n6_only,'Exclude',longidx);
n6_only_LME_names = n6_only_LME.CoefficientNames';
n6_only_LME_coefs = n6_only_LME.Coefficients.Estimate;
n6_only_LME_SE = n6_only_LME.Coefficients.SE;
n6_only_LME_Pval = n6_only_LME.Coefficients.pValue;

%% zone indicator in simple beh model
modelspec_n6_zone = 'n_duration ~ n_minus_one_Durations_All*zone_indicator+ n_minus_two_Durations_All + n_minus_three_Durations_All+ n_minus_four_Durations_All+ n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator) + (1|day_indicator)';
n6_zone_LME = fitlme(Regression_table, modelspec_n6_zone,'Exclude',longidx);


%% lmes for stim indicator
%% full model
modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints = 'n_duration ~ipi2*n_minus_two_Durations_All +   n_minus_one_Durations_All*criteria_percent_indicator + n_minus_one_Durations_All*lp_on_times + n_minus_one_Durations_All*ipi1 + Reward_idx_n1*n_minus_one_Durations_All+ HE_n1_Indicator*n_minus_one_Durations_All + n_minus_one_Durations_All * n_Stim_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  +criteria_percent_indicator + lp_on_times +moving_average_lp_length_n7andback*n_Stim_All +  moving_average_lp_length_n7andback*ipi1 + moving_average_lp_length_n7andback*HE_n1_Indicator + Reward_idx_n1*moving_average_lp_length_n7andback + moving_average_lp_length_n7andback*criteria_percent_indicator+moving_average_lp_length_n7andback*lp_on_times+(1|mouse_indicator) + (1|day_indicator)';
n6_and_stim_me_ma_ctl_zonestim_no_stimints = fitlme(Regression_table, modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints,'Exclude',longidx);
n6_and_stim_me_ma_ctl_zonestim_no_stimints_coefs = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.Estimate;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_SE = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.SE;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_pval = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.pValue;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_name = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.Name;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_df = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.DF;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_tstat = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.tStat;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_upper = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.Upper;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_lower = n6_and_stim_me_ma_ctl_zonestim_no_stimints.Coefficients.Lower;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_anova = anova(n6_and_stim_me_ma_ctl_zonestim_no_stimints);
n6_and_stim_me_ma_ctl_zonestim_no_stimints_fstat = n6_and_stim_me_ma_ctl_zonestim_no_stimints_anova.FStat;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_fpval = n6_and_stim_me_ma_ctl_zonestim_no_stimints_anova.pValue;

%HE int?
modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_he = 'n_duration ~HE_n1_Indicator:n_Stim_All +n_minus_one_Durations_All:HE_n1_Indicator:n_Stim_All + ipi2*n_minus_two_Durations_All +   n_minus_one_Durations_All*criteria_percent_indicator + n_minus_one_Durations_All*lp_on_times + n_minus_one_Durations_All*ipi1 + Reward_idx_n1*n_minus_one_Durations_All+ HE_n1_Indicator*n_minus_one_Durations_All + n_minus_one_Durations_All * n_Stim_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  +criteria_percent_indicator + lp_on_times +moving_average_lp_length_n7andback*n_Stim_All +  moving_average_lp_length_n7andback*ipi1 + moving_average_lp_length_n7andback*HE_n1_Indicator + Reward_idx_n1*moving_average_lp_length_n7andback + moving_average_lp_length_n7andback*criteria_percent_indicator+moving_average_lp_length_n7andback*lp_on_times+(1|mouse_indicator) + (1|day_indicator)';
n6_and_stim_me_ma_ctl_zonestim_no_stimints_he = fitlme(Regression_table, modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_he,'Exclude',longidx);

%ipi int?
modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_ipi = 'n_duration ~ipi1:n_Stim_All +n_minus_one_Durations_All:ipi1:n_Stim_All + ipi2*n_minus_two_Durations_All +   n_minus_one_Durations_All*criteria_percent_indicator + n_minus_one_Durations_All*lp_on_times + n_minus_one_Durations_All*ipi1 + Reward_idx_n1*n_minus_one_Durations_All+ HE_n1_Indicator*n_minus_one_Durations_All + n_minus_one_Durations_All * n_Stim_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  +criteria_percent_indicator + lp_on_times +moving_average_lp_length_n7andback*n_Stim_All +  moving_average_lp_length_n7andback*ipi1 + moving_average_lp_length_n7andback*HE_n1_Indicator + Reward_idx_n1*moving_average_lp_length_n7andback + moving_average_lp_length_n7andback*criteria_percent_indicator+moving_average_lp_length_n7andback*lp_on_times+(1|mouse_indicator) + (1|day_indicator)';
n6_and_stim_me_ma_ctl_zonestim_no_stimints_ipi= fitlme(Regression_table, modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_ipi,'Exclude',longidx);

%add in latency covariate
modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_latency = 'n_duration ~stim_and_zone_latency*n_minus_one_Durations_All*n_Stim_All+ stim_and_zone_latency*n_minus_one_Durations_All + stim_and_zone_latency*n_Stim_All +ipi2*n_minus_two_Durations_All +   n_minus_one_Durations_All*criteria_percent_indicator + n_minus_one_Durations_All*lp_on_times + n_minus_one_Durations_All*ipi1 + Reward_idx_n1*n_minus_one_Durations_All+ HE_n1_Indicator*n_minus_one_Durations_All + n_minus_one_Durations_All * n_Stim_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  +criteria_percent_indicator + lp_on_times +moving_average_lp_length_n7andback*n_Stim_All +  moving_average_lp_length_n7andback*ipi1 + moving_average_lp_length_n7andback*HE_n1_Indicator + Reward_idx_n1*moving_average_lp_length_n7andback + moving_average_lp_length_n7andback*criteria_percent_indicator+moving_average_lp_length_n7andback*lp_on_times+(1|mouse_indicator) + (1|day_indicator)';
n6_and_stim_me_ma_ctl_zonestim_no_stimints_latency = fitlme(Regression_table, modelspec_n6_and_stim_me_ma_ctl_zonestim_no_simtints_latency,'Exclude',longidx);
n6_and_stim_me_ma_ctl_zonestim_no_stimints_latency_coefs = n6_and_stim_me_ma_ctl_zonestim_no_stimints_latency.Coefficients.Estimate;
n6_and_stim_me_ma_ctl_zonestim_no_stimints_latency_SE = n6_and_stim_me_ma_ctl_zonestim_no_stimints_latency.Coefficients.SE;


%% use stim indicator to section out variables from Regression_table
Regression_table.n_Stim_All =double(Regression_table.n_Stim_All);
Regression_table.n_minus_one_Stim_All =double(Regression_table.n_minus_one_Stim_All);

stim_n_idx = Regression_table.n_Stim_All ==2;
no_stim_n_idx = Regression_table.n_Stim_All ==1;
stim_n1_idx = Regression_table.n_minus_one_Stim_All ==2;
no_stim_n1_idx = Regression_table.n_minus_one_Stim_All ==1;

stim_n_duration = Regression_table.n_duration(stim_n_idx) ;
no_stim_n_duration = Regression_table.n_duration(no_stim_n_idx) ;
%exclude anything over 10s
stim_n_duration = stim_n_duration(stim_n_duration<10000);
no_stim_n_duration = no_stim_n_duration(no_stim_n_duration<10000);

stim_n1_duration = Regression_table.n_duration(stim_n1_idx) ;
no_stim_n1_duration = Regression_table.n_duration(no_stim_n1_idx);

stim_n_minus1_duration = Regression_table.n_minus_one_Durations_All(stim_n_idx) ;
no_stim_n_minus1_duration = Regression_table.n_minus_one_Durations_All(no_stim_n_idx);
stim_n_minus1_duration = stim_n_minus1_duration(stim_n_minus1_duration<10000);
no_stim_n_minus1_duration = no_stim_n_minus1_duration(no_stim_n_minus1_duration<10000);

stim_n_ipi = Regression_table.ipi1(stim_n_idx) ;
no_stim_n_ipi = Regression_table.ipi1(no_stim_n_idx) ;

end

