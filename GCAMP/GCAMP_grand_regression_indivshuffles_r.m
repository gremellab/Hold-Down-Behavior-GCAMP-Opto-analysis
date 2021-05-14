function [] = GCAMP_grand_regression_indivshuffles_r(Grouped_GCAMP)
%Last updated 5/10/21 Drew Schreiner
%This function takes the Grouped_GCAMP data structure
%and specifically the individual lme table table
%and plops all the individual animal data together into one big table
%which is used to build the simple and complex LMEs predicting activity at
%different time points

tic
Grand_Regression_Cell =[];
Shuf_distributions = [];
Shuf_n_back_Lengths_distribution ={};

%lump all the mice together into one big cell
for i = 1:length(Grouped_GCAMP.Mice)
    Grand_Regression_Cell = [Grand_Regression_Cell; Grouped_GCAMP.Mice{1,i}.regression_cell];
    Shuf_distributions = [Shuf_distributions; Grouped_GCAMP.Mice{1,i}.Shuffled_Durs_Distribtuion ]; 
    Shuf_n_back_Lengths_distribution = [Shuf_n_back_Lengths_distribution Grouped_GCAMP.Mice{1,i}.shuf_n_back_Lengths_distribution];      
end

%get column names from the individual mouse data, since the vars dont
%change
column_names = Grouped_GCAMP.Mice{1,1}.data_table_variables;
Grand_Regression_Table = cell2table(Grand_Regression_Cell, 'VariableNames',column_names);

%add in an n zero reward, for use with post offset graphs
 Grand_Regression_Table.n_zero_reward = Grand_Regression_Table.n_minus_zero_Durations_All>=Grand_Regression_Table.criteria_indicator;

%% shuffled lp order regressions
%create vars for shuffled coefficients

shuf_lme_interp_AUC_lponly_n0_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n1_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n2_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n3_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n4_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n5_coef_ALL = [];
shuf_lme_interp_AUC_lponly_n6_coef_ALL = [];

shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL = [];
shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL = [];

shuf_lme_post_offset_base_raw_1s_n6_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n5_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n4_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n3_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n2_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n1_coef_ALL = [];
shuf_lme_post_offset_base_raw_1s_n0_coef_ALL = [];

%Find Presses over 10s, use this idx to exclude them from the model
longidx =find(Grand_Regression_Table.n_minus_zero_Durations_All>=10000);
%should be redundant, removed long presses earlier in pipeline


for j =  1:size(Shuf_n_back_Lengths_distribution,1)
%need to erase the i_shuffled array every time, this will be used to create
%the table for regression
i_shuffled =[];
for i = 1:size(Shuf_n_back_Lengths_distribution,2)
%first shuffle for all mice
i_shuffled = [i_shuffled; Shuf_n_back_Lengths_distribution{j,i}'];
end

%add in the appropriate shuffled n lp (ie, the shuffled dist that built the
%n-back array)
  i_shuffled = [i_shuffled Shuf_distributions(:,j)];
  
   T_shuffled = array2table(i_shuffled,'VariableNames',{'shuf_n_minus_one_Durations_All',...
    'shuf_n_minus_two_Durations_All', 'shuf_n_minus_three_Durations_All', 'shuf_n_minus_four_Durations_All',...
    'shuf_n_minus_five_Durations_All', 'shuf_n_minus_six_Durations_All', 'shuf_n_minus_seven_Durations_All',...
    'shuf_n_minus_eight_Durations_All', 'shuf_n_minus_nine_Durations_All', 'shuf_n_minus_ten_Durations_All',...
  'shuf_LP_Durations_All'});

T_Both = [T_shuffled Grand_Regression_Table];

%% shuffle each individual lag in the interpolated models
%e.g., only n-0 duration is shuffled
shuf_modelspec_interp_AUC_lponly_n0 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + shuf_LP_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n0 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n0,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n0_coef = shuf_lme_interp_AUC_lponly_n0.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n0_coef_ALL = [shuf_lme_interp_AUC_lponly_n0_coef_ALL shuf_lme_interp_AUC_lponly_n0_coef];

shuf_modelspec_interp_AUC_lponly_n1 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + shuf_n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n1 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n1,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n1_coef = shuf_lme_interp_AUC_lponly_n1.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n1_coef_ALL = [shuf_lme_interp_AUC_lponly_n1_coef_ALL shuf_lme_interp_AUC_lponly_n1_coef];

shuf_modelspec_interp_AUC_lponly_n2 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + shuf_n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n2 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n2,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n2_coef = shuf_lme_interp_AUC_lponly_n2.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n2_coef_ALL = [shuf_lme_interp_AUC_lponly_n2_coef_ALL shuf_lme_interp_AUC_lponly_n2_coef];

shuf_modelspec_interp_AUC_lponly_n3 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + shuf_n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n3 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n3,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n3_coef = shuf_lme_interp_AUC_lponly_n3.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n3_coef_ALL = [shuf_lme_interp_AUC_lponly_n3_coef_ALL shuf_lme_interp_AUC_lponly_n3_coef];

shuf_modelspec_interp_AUC_lponly_n4 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + shuf_n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n4 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n4,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n4_coef = shuf_lme_interp_AUC_lponly_n4.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n4_coef_ALL = [shuf_lme_interp_AUC_lponly_n4_coef_ALL shuf_lme_interp_AUC_lponly_n4_coef];

shuf_modelspec_interp_AUC_lponly_n5 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + shuf_n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n5 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n5,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n5_coef = shuf_lme_interp_AUC_lponly_n5.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n5_coef_ALL = [shuf_lme_interp_AUC_lponly_n5_coef_ALL shuf_lme_interp_AUC_lponly_n5_coef];

shuf_modelspec_interp_AUC_lponly_n6 = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + shuf_n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_interp_AUC_lponly_n6 = fitlme(T_Both,shuf_modelspec_interp_AUC_lponly_n6,'Exclude',longidx);
shuf_lme_interp_AUC_lponly_n6_coef = shuf_lme_interp_AUC_lponly_n6.Coefficients.Estimate;
shuf_lme_interp_AUC_lponly_n6_coef_ALL = [shuf_lme_interp_AUC_lponly_n6_coef_ALL shuf_lme_interp_AUC_lponly_n6_coef];

%% shuffle each lag in the pre onset
shuf_modelspec_pre_onset_base_raw_1s_n0  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + shuf_LP_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n0 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n0,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n0_coef = shuf_lme_pre_onset_base_raw_1s_n0.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL shuf_lme_pre_onset_base_raw_1s_n0_coef];

shuf_modelspec_pre_onset_base_raw_1s_n1  = ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + shuf_n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n1 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n1,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n1_coef = shuf_lme_pre_onset_base_raw_1s_n1.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL shuf_lme_pre_onset_base_raw_1s_n1_coef];

shuf_modelspec_pre_onset_base_raw_1s_n2  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + shuf_n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n2 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n2,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n2_coef = shuf_lme_pre_onset_base_raw_1s_n2.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL shuf_lme_pre_onset_base_raw_1s_n2_coef];

shuf_modelspec_pre_onset_base_raw_1s_n3  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + shuf_n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n3 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n3,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n3_coef = shuf_lme_pre_onset_base_raw_1s_n3.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL shuf_lme_pre_onset_base_raw_1s_n3_coef];

shuf_modelspec_pre_onset_base_raw_1s_n4  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + shuf_n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n4 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n4,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n4_coef = shuf_lme_pre_onset_base_raw_1s_n4.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL shuf_lme_pre_onset_base_raw_1s_n4_coef];

shuf_modelspec_pre_onset_base_raw_1s_n5  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + shuf_n_minus_five_Durations_All + n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n5 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n5,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n5_coef = shuf_lme_pre_onset_base_raw_1s_n5.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL shuf_lme_pre_onset_base_raw_1s_n5_coef];

shuf_modelspec_pre_onset_base_raw_1s_n6  =  ' base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1 + base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 +  base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + shuf_n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_pre_onset_base_raw_1s_n6 = fitlme(T_Both,shuf_modelspec_pre_onset_base_raw_1s_n6,'Exclude',longidx);
shuf_lme_pre_onset_base_raw_1s_n6_coef = shuf_lme_pre_onset_base_raw_1s_n6.Coefficients.Estimate;
shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL = [shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL shuf_lme_pre_onset_base_raw_1s_n6_coef];

%% post offset 0 to 1s
shuf_modelspec_post_offset_base_raw_1s_n0 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + shuf_LP_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n0 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n0,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n0_coef = shuf_lme_post_offset_base_raw_1s_n0.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n0_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n0_coef_ALL shuf_lme_post_offset_base_raw_1s_n0_coef];

shuf_modelspec_post_offset_base_raw_1s_n1 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + shuf_n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n1 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n1,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n1_coef = shuf_lme_post_offset_base_raw_1s_n1.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n1_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n1_coef_ALL shuf_lme_post_offset_base_raw_1s_n1_coef];

shuf_modelspec_post_offset_base_raw_1s_n2 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + shuf_n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n2 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n2,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n2_coef = shuf_lme_post_offset_base_raw_1s_n2.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n2_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n2_coef_ALL shuf_lme_post_offset_base_raw_1s_n2_coef];

shuf_modelspec_post_offset_base_raw_1s_n3 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + shuf_n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n3 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n3,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n3_coef = shuf_lme_post_offset_base_raw_1s_n3.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n3_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n3_coef_ALL shuf_lme_post_offset_base_raw_1s_n3_coef];

shuf_modelspec_post_offset_base_raw_1s_n4 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + shuf_n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n4 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n4,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n4_coef = shuf_lme_post_offset_base_raw_1s_n4.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n4_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n4_coef_ALL shuf_lme_post_offset_base_raw_1s_n4_coef];

shuf_modelspec_post_offset_base_raw_1s_n5 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + shuf_n_minus_five_Durations_All + n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n5 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n5,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n5_coef = shuf_lme_post_offset_base_raw_1s_n5.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n5_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n5_coef_ALL shuf_lme_post_offset_base_raw_1s_n5_coef];

shuf_modelspec_post_offset_base_raw_1s_n6 =  ' base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + shuf_n_minus_six_Durations_All  + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
shuf_lme_post_offset_base_raw_1s_n6 = fitlme(T_Both,shuf_modelspec_post_offset_base_raw_1s_n6,'Exclude',longidx);
shuf_lme_post_offset_base_raw_1s_n6_coef = shuf_lme_post_offset_base_raw_1s_n6.Coefficients.Estimate;
shuf_lme_post_offset_base_raw_1s_n6_coef_ALL = [shuf_lme_post_offset_base_raw_1s_n6_coef_ALL shuf_lme_post_offset_base_raw_1s_n6_coef];

end

%% first create a purely beahvioral LME predicting N duration
%n - 1 to n - 10 onlt
modelspec_behavior_lponly = 'n_minus_zero_Durations_All ~ n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All + n_minus_six_Durations_All + n_minus_seven_Durations_All + n_minus_eight_Durations_All + n_minus_nine_Durations_All + n_minus_ten_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_behavior_lponly = fitlme(Grand_Regression_Table,modelspec_behavior_lponly,'Exclude',longidx);
lme_behavior_lponly_coef = lme_behavior_lponly.Coefficients.Estimate;
lme_behavior_lponly_se = lme_behavior_lponly.Coefficients.SE;

%reduced based on model selection from 33-3
modelspec_behavior_reduced_333 = 'n_minus_zero_Durations_All ~criteria_percent_indicator+n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2+n_minus_two_Durations_All*ipi2+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+moving_average_lp_length_n7andback:criteria_percent_indicator  + lp_on_times:moving_average_lp_length_n7andback+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_behavior_reduced_333 = fitlme(Grand_Regression_Table,modelspec_behavior_reduced_333,'Exclude',longidx);
lme_behavior_reduced_coef = lme_behavior_reduced_333.Coefficients.Estimate;
lme_behavior_reduced_se = lme_behavior_reduced_333.Coefficients.SE;

%% interpolated activity during the press

%simple relation to durations
modelspec_interp_all_n6_act = 'interp_all_AUC ~ interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_interp_all_n6_act = fitlme(Grand_Regression_Table,modelspec_interp_all_n6_act,'Exclude',longidx);
lme_interp_all_n6_act_coef = lme_interp_all_n6_act.Coefficients.Estimate;
lme_interp_all_n6_act_se = lme_interp_all_n6_act.Coefficients.SE;
lme_interp_all_n6_act_pval = lme_interp_all_n6_act.Coefficients.pValue;
lme_interp_all_n6_act_name = lme_interp_all_n6_act.Coefficients.Name;
lme_interp_all_n6_act_anova = anova(lme_interp_all_n6_act);
lme_interp_all_n6_act_fstat = lme_interp_all_n6_act_anova.FStat;
lme_interp_all_n6_act_fpval = lme_interp_all_n6_act_anova.pValue;
lme_interp_all_n6_act_upper = lme_interp_all_n6_act.Coefficients.Upper;
lme_interp_all_n6_act_lower = lme_interp_all_n6_act.Coefficients.Lower;

avg_interp_AUC = nanmean(Grand_Regression_Table.interp_all_AUC);
med_interp_AUC = nanmedian(Grand_Regression_Table.interp_all_AUC);

%perm test comparing to shuffled coefs
n_p_interp_n0 = sum( abs(shuf_lme_interp_AUC_lponly_n0_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(2)))/size(shuf_lme_interp_AUC_lponly_n0_coef_ALL,2);
n_p_interp_n1 = sum( abs(shuf_lme_interp_AUC_lponly_n1_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(3)))/size(shuf_lme_interp_AUC_lponly_n1_coef_ALL,2);
n_p_interp_n2 = sum( abs(shuf_lme_interp_AUC_lponly_n2_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(4)))/size(shuf_lme_interp_AUC_lponly_n2_coef_ALL,2);
n_p_interp_n3 = sum( abs(shuf_lme_interp_AUC_lponly_n3_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(5)))/size(shuf_lme_interp_AUC_lponly_n3_coef_ALL,2);
n_p_interp_n4 = sum( abs(shuf_lme_interp_AUC_lponly_n4_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(6)))/size(shuf_lme_interp_AUC_lponly_n4_coef_ALL,2);
n_p_interp_n5 = sum( abs(shuf_lme_interp_AUC_lponly_n5_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(7)))/size(shuf_lme_interp_AUC_lponly_n5_coef_ALL,2);
n_p_interp_n6 = sum( abs(shuf_lme_interp_AUC_lponly_n6_coef_ALL(2,:)) > abs(lme_interp_all_n6_act_coef(8)))/size(shuf_lme_interp_AUC_lponly_n6_coef_ALL,2);

n_p_interp_n0_to_n6 = [n_p_interp_n0; n_p_interp_n1; n_p_interp_n2; n_p_interp_n3; n_p_interp_n4; n_p_interp_n5; n_p_interp_n6;];

%shuffled mean/sem
interp_n0_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n0_coef_ALL(2,:));
interp_n0_shuf_sem = std(shuf_lme_interp_AUC_lponly_n0_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n0_coef_ALL));
interp_n1_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n1_coef_ALL(2,:));
interp_n1_shuf_sem = std(shuf_lme_interp_AUC_lponly_n1_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n1_coef_ALL));
interp_n2_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n2_coef_ALL(2,:));
interp_n2_shuf_sem = std(shuf_lme_interp_AUC_lponly_n2_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n2_coef_ALL));
interp_n3_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n3_coef_ALL(2,:));
interp_n3_shuf_sem = std(shuf_lme_interp_AUC_lponly_n3_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n3_coef_ALL));
interp_n4_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n4_coef_ALL(2,:));
interp_n4_shuf_sem = std(shuf_lme_interp_AUC_lponly_n4_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n4_coef_ALL));
interp_n5_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n5_coef_ALL(2,:));
interp_n5_shuf_sem = std(shuf_lme_interp_AUC_lponly_n5_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n5_coef_ALL));
interp_n6_shuf_mean = mean(shuf_lme_interp_AUC_lponly_n6_coef_ALL(2,:));
interp_n6_shuf_sem = std(shuf_lme_interp_AUC_lponly_n6_coef_ALL(2,:))/sqrt(length(shuf_lme_interp_AUC_lponly_n6_coef_ALL));

interp_n0_n6_mean = [interp_n0_shuf_mean; interp_n1_shuf_mean; interp_n2_shuf_mean; interp_n3_shuf_mean; interp_n4_shuf_mean; interp_n5_shuf_mean; interp_n6_shuf_mean ];
interp_n0_n6_sem = [interp_n0_shuf_sem; interp_n1_shuf_sem; interp_n2_shuf_sem; interp_n3_shuf_sem; interp_n4_shuf_sem; interp_n5_shuf_sem; interp_n6_shuf_sem ];

%now add in the cannonical beh model terms in addition
modelspec_interp_all_n6_act_andints = 'interp_all_AUC ~interp_all_AUC_n1 + interp_all_AUC_n2 + interp_all_AUC_n3 + interp_all_AUC_n4 + interp_all_AUC_n5 + interp_all_AUC_n6 + n_minus_zero_Durations_All + criteria_percent_indicator+n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+moving_average_lp_length_n7andback:criteria_percent_indicator  + lp_on_times:moving_average_lp_length_n7andback+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_interp_all_n6_act_andints = fitlme(Grand_Regression_Table,modelspec_interp_all_n6_act_andints,'Exclude',longidx)
lme_interp_all_n6_act_andints_names = lme_interp_all_n6_act_andints.Coefficients.Name;
lme_interp_all_n6_act_andints_coef = lme_interp_all_n6_act_andints.Coefficients.Estimate;
lme_interp_all_n6_act_andints_pvals = lme_interp_all_n6_act_andints.Coefficients.pValue;	
lme_interp_all_n6_act_andints_se = lme_interp_all_n6_act_andints.Coefficients.SE;	
lme_interp_all_n6_act_andints_df = lme_interp_all_n6_act_andints.Coefficients.DF;	
lme_interp_all_n6_act_andints_upper = lme_interp_all_n6_act_andints.Coefficients.Upper;	
lme_interp_all_n6_act_andints_lower = lme_interp_all_n6_act_andints.Coefficients.Lower;	
lme_interp_all_n6_act_andints_anova = anova(lme_interp_all_n6_act_andints);
lme_interp_all_n6_act_andints_anova_Fval = lme_interp_all_n6_act_andints_anova.FStat;
lme_interp_all_n6_act_andints_anova_Fpval = lme_interp_all_n6_act_andints_anova.pValue;

%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_interp_all_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.interp_all_AUC <=yhatCI(:,2) &  Grand_Regression_Table.interp_all_AUC >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.interp_all_AUC)

 %M2
 %23.7
 %M2-DMS
 %21.7
 
[yhat yhatCI] =predict(lme_interp_all_n6_act_andints,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.interp_all_AUC <=yhatCI(:,2) &  Grand_Regression_Table.interp_all_AUC >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.interp_all_AUC)

 %M2
 %32.7%, 9% better than simple model
 %M2-DMS
 %31.2%, 9.5% better than simple model

%% pre onset 1s 
modelspec_pre_onset_base_raw_1s_n6_act = 'base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1+ base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 + base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All + (1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_pre_onset_base_raw_1s_n6_act = fitlme(Grand_Regression_Table,modelspec_pre_onset_base_raw_1s_n6_act,'Exclude',longidx);
lme_pre_onset_base_raw_1s_n6_act_coef = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Estimate;
lme_pre_onset_base_raw_1s_n6_act_se = lme_pre_onset_base_raw_1s_n6_act.Coefficients.SE;
lme_pre_onset_base_raw_1s_n6_act_pval = lme_pre_onset_base_raw_1s_n6_act.Coefficients.pValue;
lme_pre_onset_base_raw_1s_n6_act_name = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Name;
lme_pre_onset_base_raw_1s_n6_act_anova = anova(lme_pre_onset_base_raw_1s_n6_act);
lme_pre_onset_base_raw_1s_n6_act_fstat = lme_pre_onset_base_raw_1s_n6_act_anova.FStat;
lme_pre_onset_base_raw_1s_n6_act_fpval = lme_pre_onset_base_raw_1s_n6_act_anova.pValue;
lme_pre_onset_base_raw_1s_n6_act_upper = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Upper;
lme_pre_onset_base_raw_1s_n6_act_lower = lme_pre_onset_base_raw_1s_n6_act.Coefficients.Lower;

avg_preonset_AUC = nanmean(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw);
med_preonset_AUC = nanmedian(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw);

n_p_pre_n0 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(2)))/size(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL,2);
n_p_pre_n1 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(3)))/size(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL,2);
n_p_pre_n2 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(4)))/size(shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL,2);
n_p_pre_n3 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(5)))/size(shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL,2);
n_p_pre_n4 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(6)))/size(shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL,2);
n_p_pre_n5 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(7)))/size(shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL,2);
n_p_pre_n6 = sum( abs(shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL(2,:)) > abs(lme_pre_onset_base_raw_1s_n6_act_coef(8)))/size(shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL,2);

n_p_pre_n0_to_n6 = [n_p_pre_n0; n_p_pre_n1; n_p_pre_n2; n_p_pre_n3; n_p_pre_n4; n_p_pre_n5; n_p_pre_n6;];

pre_n0_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL(2,:));
pre_n0_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n0_coef_ALL));
pre_n1_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL(2,:));
pre_n1_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n1_coef_ALL));
pre_n2_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL(2,:));
pre_n2_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n2_coef_ALL));
pre_n3_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL(2,:));
pre_n3_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n3_coef_ALL));
pre_n4_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL(2,:));
pre_n4_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n4_coef_ALL));
pre_n5_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL(2,:));
pre_n5_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n5_coef_ALL));
pre_n6_shuf_mean = mean(shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL(2,:));
pre_n6_shuf_sem = std(shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL(2,:))/sqrt(length(shuf_lme_pre_onset_base_raw_1s_n6_coef_ALL));

pre_n0_n6_mean = [pre_n0_shuf_mean; pre_n1_shuf_mean; pre_n2_shuf_mean; pre_n3_shuf_mean; pre_n4_shuf_mean; pre_n5_shuf_mean; pre_n6_shuf_mean ];
pre_n0_n6_sem = [pre_n0_shuf_sem; pre_n1_shuf_sem; pre_n2_shuf_sem; pre_n3_shuf_sem; pre_n4_shuf_sem; pre_n5_shuf_sem; pre_n6_shuf_sem ];

%pre onset 1s n1n6 ints
modelspec_preonset_all_n6_act_andints = 'base_neg_1_to_onset_AUC_raw ~ base_neg_1_to_onset_AUC_raw_n1+ base_neg_1_to_onset_AUC_raw_n2 + base_neg_1_to_onset_AUC_raw_n3 + base_neg_1_to_onset_AUC_raw_n4 + base_neg_1_to_onset_AUC_raw_n5 + base_neg_1_to_onset_AUC_raw_n6 + n_minus_zero_Durations_All + criteria_percent_indicator+n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+moving_average_lp_length_n7andback:criteria_percent_indicator  + lp_on_times:moving_average_lp_length_n7andback+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_preonset_all_n6_act_andints = fitlme(Grand_Regression_Table,modelspec_preonset_all_n6_act_andints,'Exclude',longidx)
lme_preonset_all_n6_act_andints_coef = lme_preonset_all_n6_act_andints.Coefficients.Estimate;
lme_preonset_all_n6_act_andints_se = lme_preonset_all_n6_act_andints.Coefficients.SE;
lme_preonset_all_n6_act_andints_pvals = lme_preonset_all_n6_act_andints.Coefficients.pValue;	
lme_preonset_all_n6_act_andints_names = lme_preonset_all_n6_act_andints.Coefficients.Name;
lme_preonset_all_n6_act_andints_df = lme_preonset_all_n6_act_andints.Coefficients.DF;	
lme_preonset_all_n6_act_andints_upper = lme_preonset_all_n6_act_andints.Coefficients.Upper;	
lme_preonset_all_n6_act_andints_lower = lme_preonset_all_n6_act_andints.Coefficients.Lower;	
lme_preonset_all_n6_act_andints_anova = anova(lme_preonset_all_n6_act_andints);
lme_preonset_all_n6_act_andints_anova_Fval = lme_preonset_all_n6_act_andints_anova.FStat;
lme_preonset_all_n6_act_andints_anova_Fpval = lme_preonset_all_n6_act_andints_anova.pValue;

%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_pre_onset_base_raw_1s_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_neg_1_to_onset_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_neg_1_to_onset_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw)

 %M2
 %36.4%
 %M2-DMS
 %28.1%
 
[yhat yhatCI] =predict(lme_preonset_all_n6_act_andints,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_neg_1_to_onset_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_neg_1_to_onset_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_neg_1_to_onset_AUC_raw)

 %M2
 %50.0%, 13.6% better than simple
 %M2-DMS
 %41.9%, 13.8% better than simple model

%% post offset 1s
modelspec_post_offset_base_raw_to1_n6_act = 'base_offset_to1_AUC_raw ~ base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_post_offset_base_raw_to1_n6_act = fitlme(Grand_Regression_Table,modelspec_post_offset_base_raw_to1_n6_act,'Exclude',longidx);
lme_post_offset_base_raw_to1_n6_act_name = lme_post_offset_base_raw_to1_n6_act.Coefficients.Name;
lme_post_offset_base_raw_to1_n6_act_coef = lme_post_offset_base_raw_to1_n6_act.Coefficients.Estimate;
lme_post_offset_base_raw_to1_n6_act_se = lme_post_offset_base_raw_to1_n6_act.Coefficients.SE;
lme_post_offset_base_raw_to1_n6_act_pval = lme_post_offset_base_raw_to1_n6_act.Coefficients.pValue;
lme_post_offset_base_raw_to1_n6_act_name = lme_post_offset_base_raw_to1_n6_act.Coefficients.Name;
lme_post_offset_base_raw_to1_n6_act_anova = anova(lme_post_offset_base_raw_to1_n6_act);
lme_post_offset_base_raw_to1_n6_act_fstat = lme_post_offset_base_raw_to1_n6_act_anova.FStat;
lme_post_offset_base_raw_to1_n6_act_fpval = lme_post_offset_base_raw_to1_n6_act_anova.pValue;
lme_post_offset_base_raw_to1_n6_act_upper = lme_post_offset_base_raw_to1_n6_act.Coefficients.Upper;
lme_post_offset_base_raw_to1_n6_act_lower = lme_post_offset_base_raw_to1_n6_act.Coefficients.Lower;

avg_post_offset_AUC = nanmean(Grand_Regression_Table.base_offset_to1_AUC_raw);
med_post_offset_AUC = nanmedian(Grand_Regression_Table.base_offset_to1_AUC_raw);

n_p_post_n0 = sum( abs(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(2)))/size(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL,2);
n_p_post_n1 = sum( abs(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(3)))/size(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL,2);
n_p_post_n2 = sum( abs(shuf_lme_post_offset_base_raw_1s_n2_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(4)))/size(shuf_lme_post_offset_base_raw_1s_n2_coef_ALL,2);
n_p_post_n3 = sum( abs(shuf_lme_post_offset_base_raw_1s_n3_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(5)))/size(shuf_lme_post_offset_base_raw_1s_n3_coef_ALL,2);
n_p_post_n4 = sum( abs(shuf_lme_post_offset_base_raw_1s_n4_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(6)))/size(shuf_lme_post_offset_base_raw_1s_n4_coef_ALL,2);
n_p_post_n5 = sum( abs(shuf_lme_post_offset_base_raw_1s_n5_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(7)))/size(shuf_lme_post_offset_base_raw_1s_n5_coef_ALL,2);
n_p_post_n6 = sum( abs(shuf_lme_post_offset_base_raw_1s_n6_coef_ALL(2,:)) > abs(lme_post_offset_base_raw_to1_n6_act_coef(8)))/size(shuf_lme_post_offset_base_raw_1s_n6_coef_ALL,2);

n_p_post_n0_to_n6 = [n_p_post_n0; n_p_post_n1; n_p_post_n2; n_p_post_n3; n_p_post_n4; n_p_post_n5; n_p_post_n6;];

post_n0_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL(2,:));
post_n0_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n0_coef_ALL));
post_n1_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL(2,:));
post_n1_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n1_coef_ALL));
post_n2_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n2_coef_ALL(2,:));
post_n2_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n2_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n2_coef_ALL));
post_n3_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n3_coef_ALL(2,:));
post_n3_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n3_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n3_coef_ALL));
post_n4_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n4_coef_ALL(2,:));
post_n4_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n4_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n4_coef_ALL));
post_n5_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n5_coef_ALL(2,:));
post_n5_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n5_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n5_coef_ALL));
post_n6_shuf_mean = mean(shuf_lme_post_offset_base_raw_1s_n6_coef_ALL(2,:));
post_n6_shuf_sem = std(shuf_lme_post_offset_base_raw_1s_n6_coef_ALL(2,:))/sqrt(length(shuf_lme_post_offset_base_raw_1s_n6_coef_ALL));

post_n0_n6_mean = [post_n0_shuf_mean; post_n1_shuf_mean; post_n2_shuf_mean; post_n3_shuf_mean; post_n4_shuf_mean; post_n5_shuf_mean; post_n6_shuf_mean ];
post_n0_n6_sem = [post_n0_shuf_sem; post_n1_shuf_sem; post_n2_shuf_sem; post_n3_shuf_sem; post_n4_shuf_sem; post_n5_shuf_sem; post_n6_shuf_sem ];

% add in n_zero_reward (i.e., was the just finished lp rewarded) and
% interaction
% n_minus_zero_Durations_All*n_zero_reward +
% moving_average_lp_length_n7andback:n_zero_reward
modelspec_postoffset_all_n6_act_andints_n0 = 'base_offset_to1_AUC_raw ~n_zero_reward:moving_average_lp_length_n7andback+ n_minus_zero_Durations_All * n_zero_reward + base_offset_to1_AUC_raw_n1 + base_offset_to1_AUC_raw_n2 + base_offset_to1_AUC_raw_n3 + base_offset_to1_AUC_raw_n4 + base_offset_to1_AUC_raw_n5 + base_offset_to1_AUC_raw_n6 + n_minus_zero_Durations_All + criteria_percent_indicator+n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+moving_average_lp_length_n7andback:criteria_percent_indicator  + lp_on_times:moving_average_lp_length_n7andback+(1|mouse_indicator)+(1|day_indicator)+(1|criteria_indicator)';
lme_postoffset_all_n6_act_andints_n0 = fitlme(Grand_Regression_Table,modelspec_postoffset_all_n6_act_andints_n0,'Exclude',longidx)
lme_postoffset_all_n6_act_andints_n0_coef = lme_postoffset_all_n6_act_andints_n0.Coefficients.Estimate;
lme_postoffset_all_n6_act_andints_n0_pvals = lme_postoffset_all_n6_act_andints_n0.Coefficients.pValue;	
lme_postoffset_all_n6_act_andints_n0_se = lme_postoffset_all_n6_act_andints_n0.Coefficients.SE;	
lme_postoffset_all_n6_act_andints_n0_names = lme_postoffset_all_n6_act_andints_n0.Coefficients.Name;
lme_postoffset_all_n6_act_andints_n0_df = lme_postoffset_all_n6_act_andints_n0.Coefficients.DF;	
lme_postoffset_all_n6_act_andints_n0_upper = lme_postoffset_all_n6_act_andints_n0.Coefficients.Upper;	
lme_postoffset_all_n6_act_andints_n0_lower = lme_postoffset_all_n6_act_andints_n0.Coefficients.Lower;	
lme_postoffset_all_n6_act_andints_n0_anova = anova(lme_postoffset_all_n6_act_andints_n0);
lme_postoffset_all_n6_act_andints_n0_anova_Fval = lme_postoffset_all_n6_act_andints_n0_anova.FStat;
lme_postoffset_all_n6_act_andints_n0_anova_Fpval = lme_postoffset_all_n6_act_andints_n0_anova.pValue;

%How well do the simple/complex models predict activity?
[yhat yhatCI] =predict(lme_post_offset_base_raw_to1_n6_act,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_offset_to1_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_offset_to1_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_offset_to1_AUC_raw)

 %M2
 %36.5%
 %M2-DMS
 %25.5%
 
[yhat yhatCI] =predict(lme_postoffset_all_n6_act_andints_n0,Grand_Regression_Table,'Simultaneous',true);
correctish_pred =Grand_Regression_Table.base_offset_to1_AUC_raw <=yhatCI(:,2) &  Grand_Regression_Table.base_offset_to1_AUC_raw >=yhatCI(:,1);
correctCI_prop = sum(correctish_pred)/length(Grand_Regression_Table.base_offset_to1_AUC_raw)

  %M2
 %53.5%, 17% better than simple
 %M2-DMS
 %41.1%, 15.6% better than simple model


%% Optional code to save the workspace
% save('1000shufdata-m2-1600all' , '-v7.3')
toc
end

