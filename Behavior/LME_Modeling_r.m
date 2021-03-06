%% Purpose
%last edited 5/10/21 Drew Schreiner
%The main purpose of this script is to create Linear Mixed Effect Models (LMEs)
%to predict the duration of lever presses (n) given subjective experiential
%information. This includes things like n-back press durations, n-back
%headentries (HE), n-back reward, interpressinterval (IPI), time in
%session, % of presses over criteria.
%These are mixed effect models becuase they will also include random
%effects of day of training and mouse, to account for the repeated nature
%of the datasets. There are sections here for creating relatively simple
%LMEs, as well as more complex LMEs that are arrived at via a priori
%hypothesis/questions of interest combined with model selection.

%an example simple model would seek to predict n duration given prior press
%durations, whether those presses were rewarded, and interactions between
%those terms
%a = regression coeffecient
%E = error term
% C =  press duration
%R = reward or no
%Cn = ao +  + a1C(n-1) + a2C(n-2) + a3C(n-3)

%           + a4R(n-1) + a5R(n-2) + a6R(n-3)
%           + a7C(n-1)R(n-1) + a8C(n-2)R(n-2) +
%            a9C(n-3)R(n-3) + E(t)
 
%R will be coded as 0 for failed press, or 1 for rewarded
%C will be duration in ms

%There is also some other code for smaller questions that will still
%utilize the Data structure obtained from the
%MEDPC_Behavior_Extract_For_Regressions_r
%This includes things like generating within session performance graphs,
%histograms of press durations, and raw duration
%values before, during, or after optogenetic stimulation (if applicable).

%For this, you need to have a Data file generated from MEDPC_Behavior_Extract_For_Regressions
%with this Data structure being loaded under "load" (you must be in the
%folder with the data structure).
%The Data File is Organized First by Day of Training: Data.Day
%The Next level is mouse: Data.Day(x).Mouse
%Within mouse, we have a large structure that includes all the variables for
%all the mice on a given day. Some of these are single numbers (e.g., total number of lever presses)
%while others are cell arrays (e.g. all the individual lever presses a
%mouse made).
%This table goes in "mouse order". So whichever mouse was first
%(numerically) in the original folder will be first here. 
%EACH ROW IS ONE MOUSE.
%The Data.Day(x).Mouse.Name variable will label each mouse.


%% Load the Data and create variables
tic
clear all
% close all

load('33-3-beh-revised-10shuf') %load the data structure

total_mice = size(Data.Day(2).Mouse,2);
%get number of days
day_number = size(Data.Day,2);
 
%Initialize Variables
LP_Durations_Cell = {};
LP_Timestamps_Cell ={};
Rewarded_Durations_Cell = {};
Unrewarded_Durations_Cell ={};
IPI_Cell ={};
Shuffled_LP_Duration_Cell ={};
HE_Indicator_Cell ={};
mouse_indicator = [];
day_indicator=[];
Logcial_LP_Cell ={};
cell_length_var = 0;
Shufffled_Logical_Lever_Press_Cell = {};
HE_n_1_Indicator_Cell ={};
HE_n_2_Indicator_Cell = {};
HE_n_3_Indicator_Cell = {};
HE_n_4_Indicator_Cell = {};
HE_n_5_Indicator_Cell = {};
HE_n_6_Indicator_Cell = {};
IPI_nback_Cell ={};
criteria_indicator_Cell={};
criteria_percent_indicator_Cell={};
moving_mean_cell ={};
after_success_duration_Cell ={};
after_failure_duration_Cell = {};
after_success_change_Cell = {};
after_failure_change_Cell = {};
up_state_idx_n1_Cell ={};
up_state_idx_n2_Cell ={};
up_state_idx_n3_Cell ={};
up_state_idx_n4_Cell ={};
up_state_idx_n5_Cell ={};
up_state_idx_n6_Cell ={};
avg_duration_indicator_Cell ={};

reward_indicator_Cell ={};
reward_indicator_n1_Cell ={};
reward_indicator_n2_Cell ={};
reward_indicator_n3_Cell ={};
reward_indicator_n4_Cell ={};
moving_average_lp_length_n7andback_Cell ={};

%used with opto stim
stim_indicator_cell = {};
stim_indicator_n1_cell = {};
stim_indicator_n2_cell = {};
stim_indicator_n3_cell = {};
stim_indicator_n4_cell = {};
stim_indicator_n5_cell = {};
stim_indicator_n6_cell = {};

histogram_counts_Cell ={};
histogram_bins_Cell={};
histogam_indexs_Cell={};
Up_State_Consecutive_Mean_All =[]; 
Up_state_prop_All=[];
Up_State_Consecutive_Mean_All_Shuf_dist_mean =[]; 
Up_state_prop_All_Shuf_dist_mean=[];
up_state_idx_cell ={};
uppersum_sem_cell ={};
Up_State_Consecutive_Mean_All_Shuf_dist =[];
Up_state_prop_All_Shuf_dist=[];

%indiv lme predictions
indiv_correctCIprop_criteria_Cell ={};
indiv_correctCIprop_criteria_complex_Cell ={};
indiv_lme_10_indiv_se_Cell ={};
indiv_lme_10_indiv_coef_Cell ={};
indiv_lme_reduced_se_Cell ={};
indiv_lme_reduced_coef_Cell ={};
%individual lme r^2s
r_squared_criteria_Cell={};
r_squared_adjusted_criteria_Cell = {};
r_squared_complex_criteria_Cell = {};
r_squared_adjusted_complex_criteria_Cell = {};
mouse_number_for_corrs =[];
day_number_for_corrs =[];
 
mean_Shuffled_Length_After_Success_All_Shuffles =[];
mean_Shuffled_Length_After_Failure_All_Shuffles =[];
std_Shuffled_Length_After_Success_All_Shuffles =[];
std_Shuffled_Length_After_Failure_All_Shuffles =[];
mean_Shuffled_Change_After_Success_All_Shuffles =[];
mean_Shuffled_Change_After_Failure_All_Shuffles =[];

%% Loop through data structure to get data input into Cells
for i = 1:day_number %outer loop goes across days  
    for j = 1:total_mice %inner loop goes across mice within a day        
    
     IPI_Cell = [IPI_Cell Data.Day(i).Mouse(j).Lever_Press_IPI];
     mean_Shuffled_Length_After_Success_All_Shuffles = [mean_Shuffled_Length_After_Success_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Length_After_Success_Dist];
     mean_Shuffled_Length_After_Failure_All_Shuffles = [mean_Shuffled_Length_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Length_After_Failure_Dist];
     std_Shuffled_Length_After_Success_All_Shuffles = [std_Shuffled_Length_After_Success_All_Shuffles; Data.Day(i).Mouse(j).std_Shuffled_Lever_Press_Length_After_Success_Dist];
     std_Shuffled_Length_After_Failure_All_Shuffles = [std_Shuffled_Length_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).std_Shuffled_Lever_Press_Length_After_Failure_Dist];
     mean_Shuffled_Change_After_Success_All_Shuffles = [mean_Shuffled_Change_After_Success_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Change_After_Success_Dist];
     mean_Shuffled_Change_After_Failure_All_Shuffles = [mean_Shuffled_Change_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Change_After_Failure_Dist]; 
     mouse_number_for_corrs = [mouse_number_for_corrs; j];
     day_number_for_corrs = [day_number_for_corrs; i];
     indiv_lme_10_indiv_se_Cell = [indiv_lme_10_indiv_se_Cell Data.Day(i).Mouse(j).lme_10_indiv_se];
     indiv_lme_10_indiv_coef_Cell = [indiv_lme_10_indiv_coef_Cell   Data.Day(i).Mouse(j).lme_10_indiv_coef];
     indiv_lme_reduced_se_Cell =[indiv_lme_reduced_se_Cell Data.Day(i).Mouse(j).lme__moderate_remove_se];
     indiv_lme_reduced_coef_Cell =[indiv_lme_reduced_coef_Cell Data.Day(i).Mouse(j).lme__moderate_remove_coef];
     indiv_correctCIprop_criteria_Cell =[indiv_correctCIprop_criteria_Cell Data.Day(i).Mouse(j).correctCIprop_criteria];
     indiv_correctCIprop_criteria_complex_Cell =[indiv_correctCIprop_criteria_complex_Cell Data.Day(i).Mouse(j).correctCIprop_criteria_complex];
     r_squared_criteria_Cell=[r_squared_criteria_Cell  Data.Day(i).Mouse(j).r_squared_criteria];
     r_squared_adjusted_criteria_Cell = [r_squared_adjusted_criteria_Cell Data.Day(i).Mouse(j).r_squared_adjusted_criteria];
     r_squared_complex_criteria_Cell = [r_squared_complex_criteria_Cell Data.Day(i).Mouse(j).r_squared_complex_criteria ];
     r_squared_adjusted_complex_criteria_Cell = [r_squared_adjusted_complex_criteria_Cell Data.Day(i).Mouse(j).r_squared_adjusted_complex_criteria ];  
     LP_Durations_Cell = [LP_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_duration];
     Logcial_LP_Cell = [Logcial_LP_Cell Data.Day(i).Mouse(j).Logical_LPs];
     LP_Timestamps_Cell = [LP_Timestamps_Cell Data.Day(i).Mouse(j).Lever_Press_ts];
     Rewarded_Durations_Cell = [Rewarded_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_Rewarded_Lengths];
     Unrewarded_Durations_Cell = [Unrewarded_Durations_Cell Data.Day(i).Mouse(j).Lever_Press_Not_Rewarded_Lengths];
     moving_mean_cell=[moving_mean_cell Data.Day(i).Mouse(j).moving_average_lp_length];
     moving_average_lp_length_n7andback_Cell = [moving_average_lp_length_n7andback_Cell Data.Day(i).Mouse(j).moving_average_lp_length_n7andback];
     %create indicator variables that will give us an array to mark every LP
     %with the mouse and day it happened. These will be the Random Effects
     cell_length_var = cell_length_var + 1;
     %we multiple the ones here by j, the loop iteration variable for mice.
     %thus mice, won't have the same numbers in the folders, but will
     %instead be labeled 1:total_mice based on their order in the Data
     %Structure. Ditto for day_indicator, where we multiple an array of
     %ones by i (the day iteration variable). Of note, these indicators are
     %made to be as long as the Lever press durations, so that EACH LP is
     %tagged both with mouse and day.
     mouse_indicator = [mouse_indicator; (ones(size(LP_Durations_Cell{:,cell_length_var}))*j)];
     day_indicator = [day_indicator;     (ones(size(LP_Durations_Cell{:,cell_length_var}))*i)];
     %these indicators already have the appropriate size
     criteria_indicator_Cell = [criteria_indicator_Cell Data.Day(i).Mouse(j).criteria_indicator];
     criteria_percent_indicator_Cell = [criteria_percent_indicator_Cell  Data.Day(i).Mouse(j).crit_percent_indicator];
     avg_duration_indicator_Cell = [avg_duration_indicator_Cell  Data.Day(i).Mouse(j).avg_duration_indicator];
     reward_indicator_Cell =[reward_indicator_Cell Data.Day(i).Mouse(j).RE_Indicator];
     reward_indicator_n1_Cell =[reward_indicator_n1_Cell Data.Day(i).Mouse(j).RE_n_1_Indicator];
     reward_indicator_n2_Cell =[reward_indicator_n2_Cell Data.Day(i).Mouse(j).RE_n_2_Indicator];
     reward_indicator_n3_Cell =[reward_indicator_n3_Cell Data.Day(i).Mouse(j).RE_n_3_Indicator];
     reward_indicator_n4_Cell =[reward_indicator_n4_Cell Data.Day(i).Mouse(j).RE_n_4_Indicator];
     up_state_idx_cell =[up_state_idx_cell Data.Day(i).Mouse(j).Up_State_idx_Logical];
     uppersum_sem_cell =[uppersum_sem_cell Data.Day(i).Mouse(j).uppersum_sem];
     
     stim_indicator_cell = [stim_indicator_cell Data.Day(i).Mouse(j).stim_ind_all];
     stim_indicator_n1_cell =[stim_indicator_n1_cell Data.Day(i).Mouse(j).stim_ind_all_n1];
     stim_indicator_n2_cell = [stim_indicator_n2_cell Data.Day(i).Mouse(j).stim_ind_all_n2 ];
     stim_indicator_n3_cell  = [stim_indicator_n3_cell Data.Day(i).Mouse(j).stim_ind_all_n3 ];
     stim_indicator_n4_cell  = [stim_indicator_n4_cell Data.Day(i).Mouse(j).stim_ind_all_n4 ];
     stim_indicator_n5_cell  = [stim_indicator_n5_cell  Data.Day(i).Mouse(j).stim_ind_all_n5 ];
     stim_indicator_n6_cell  = [stim_indicator_n6_cell Data.Day(i).Mouse(j).stim_ind_all_n6 ];
     
     Shuffled_LP_Duration_Cell =[Shuffled_LP_Duration_Cell Data.Day(i).Mouse(j).Shuffled_Lever_Press_lengths];
     HE_Indicator_Cell = [HE_Indicator_Cell Data.Day(i).Mouse(j).HE_Indicator];
     Shufffled_Logical_Lever_Press_Cell =[Shufffled_Logical_Lever_Press_Cell Data.Day(i).Mouse(j).Shuffled_Logical_Lever_press];
     HE_n_1_Indicator_Cell = [HE_n_1_Indicator_Cell Data.Day(i).Mouse(j).HE_n_1_Indicator];
     HE_n_2_Indicator_Cell = [HE_n_2_Indicator_Cell Data.Day(i).Mouse(j).HE_n_2_Indicator];
     HE_n_3_Indicator_Cell = [HE_n_3_Indicator_Cell Data.Day(i).Mouse(j).HE_n_3_Indicator];
     HE_n_4_Indicator_Cell = [HE_n_4_Indicator_Cell Data.Day(i).Mouse(j).HE_n_4_Indicator];
     HE_n_5_Indicator_Cell = [HE_n_5_Indicator_Cell Data.Day(i).Mouse(j).HE_n_5_Indicator];
     HE_n_6_Indicator_Cell = [HE_n_6_Indicator_Cell Data.Day(i).Mouse(j).HE_n_6_Indicator];
     after_success_duration_Cell = [after_success_duration_Cell Data.Day(i).Mouse(j).Lever_Press_Length_After_Success];
     after_failure_duration_Cell = [after_failure_duration_Cell Data.Day(i).Mouse(j).Lever_Press_Length_After_Failure];
     after_success_change_Cell = [after_success_change_Cell Data.Day(i).Mouse(j).Lever_Press_Length_After_Success_Change];
     after_failure_change_Cell = [after_failure_change_Cell Data.Day(i).Mouse(j).Lever_Press_Length_After_Failure_Change];
     up_state_idx_n1_Cell =[up_state_idx_n1_Cell Data.Day(i).Mouse(j).Up_State_idx_n1];
     up_state_idx_n2_Cell =[up_state_idx_n2_Cell Data.Day(i).Mouse(j).Up_State_idx_n2];
     up_state_idx_n3_Cell =[up_state_idx_n3_Cell Data.Day(i).Mouse(j).Up_State_idx_n3];
     up_state_idx_n4_Cell =[up_state_idx_n4_Cell Data.Day(i).Mouse(j).Up_State_idx_n4];
     up_state_idx_n5_Cell =[up_state_idx_n5_Cell Data.Day(i).Mouse(j).Up_State_idx_n5];
     up_state_idx_n6_Cell =[up_state_idx_n6_Cell Data.Day(i).Mouse(j).Up_State_idx_n6]; 
     histogram_counts_Cell =[histogram_counts_Cell  Data.Day(i).Mouse(j).histogram_counts];
     histogram_bins_Cell=[histogram_bins_Cell Data.Day(i).Mouse(j).histogram_edges];
     histogam_indexs_Cell=[histogam_indexs_Cell Data.Day(i).Mouse(j).histogram_bin_indexs];
     Up_State_Consecutive_Mean_All =[Up_State_Consecutive_Mean_All  Data.Day(i).Mouse(j).Up_State_Consecutive_mean]; 
     Up_state_prop_All=[Up_state_prop_All Data.Day(i).Mouse(j).up_state_prop];  
     Up_State_Consecutive_Mean_All_Shuf_dist_mean =[Up_State_Consecutive_Mean_All_Shuf_dist_mean Data.Day(i).Mouse(j).Up_State_Consecutive_Shuf_Mean_Dist_Mean ]; 
     Up_state_prop_All_Shuf_dist_mean=[Up_state_prop_All_Shuf_dist_mean Data.Day(i).Mouse(j).Up_State_Prop_Above_Shuf_Dist_Mean];
     Up_State_Consecutive_Mean_All_Shuf_dist =[Up_State_Consecutive_Mean_All_Shuf_dist; Data.Day(i).Mouse(j).Up_State_Consecutive_Shuf_Mean ]; 
     Up_state_prop_All_Shuf_dist=[Up_state_prop_All_Shuf_dist; Data.Day(i).Mouse(j).Up_State_Prop_Above_Shuf];
          
    end %mouse loop end   
    
end %day loop end
 
%% Look at durations and statistics based on earning/not earning a reward
%Part of the purpose here is to look for any mice that develop stereotypies
%Perhaps the best indicator is the comparison of SD after a successful
%press in comparison to shuffled data

%reshape the cells to get them into cells that are num mice(row) by num day
%(col). This is used primarilyfor graphing purposes
after_failure_duration_By_Mouse = reshape(after_failure_duration_Cell,[total_mice,day_number]);
after_success_duration_By_Mouse = reshape(after_success_duration_Cell,[total_mice,day_number]);
after_failure_duration_mean = cellfun(@nanmean,after_failure_duration_By_Mouse);
after_failure_duration_std = cellfun(@nanstd,after_failure_duration_By_Mouse);
after_failure_duration_std_overallmouse = mean(after_failure_duration_std,2);
after_success_duration_mean =  cellfun(@nanmean,after_success_duration_By_Mouse);
after_success_duration_mean_overallmouse = mean(after_success_duration_mean,2);
after_success_duration_mean_grandmean = mean(after_success_duration_mean_overallmouse([1:4 6:end]));
after_success_duration_std =  cellfun(@nanstd,after_success_duration_By_Mouse);
after_success_duration_std_overallmouse = mean(after_success_duration_std,2);
after_success_duration_std_grandavg = mean(after_success_duration_std_overallmouse([1:4 6:end]));

after_success_duration_std_overallmouse_stdofdays =std(after_success_duration_std,0,2);

success_failure_std_diff = after_success_duration_std - after_failure_duration_std;
success_failure_std_diff_overallmouse = nanmean(success_failure_std_diff,2);
success_failurestd_diff_negative =success_failure_std_diff < 0;
prop_neg_success_failure_std_diff = sum(success_failurestd_diff_negative,2)/length(success_failurestd_diff_negative);

success_failure_mean_diff = after_success_duration_mean-after_failure_duration_mean;
success_failure_mean_diff_overallmouse = nanmean(success_failure_mean_diff,2);
success_failure_std_diff_overallmouse = std(success_failure_mean_diff,0,2);
%calculaute number of days that were negative per mouse
success_failure_mean_diff_negative =success_failure_mean_diff < 0;
prop_neg_success_failure_diff = sum(success_failure_mean_diff_negative,2)/length(success_failure_mean_diff_negative);

%shuffled after success/failure
overall_shuffled_after_success_failure_diff = mean_Shuffled_Length_After_Success_All_Shuffles-mean_Shuffled_Length_After_Failure_All_Shuffles;
overall_shuffled_after_success_mean = mean(mean_Shuffled_Length_After_Success_All_Shuffles,2);
 average_shuffled_after_success_std = mean(std_Shuffled_Length_After_Success_All_Shuffles,2); 
average_shuffled_after_success_std_byday = reshape(average_shuffled_after_success_std ,[total_mice,day_number]);
average_shuffled_after_success_std_bymouse =mean(average_shuffled_after_success_std_byday,2);
average_shuffled_after_success_std_bymouse_stdofdays =std(average_shuffled_after_success_std_byday,0,2);
after_success_std_actual_shuffle_diff = after_success_duration_std - average_shuffled_after_success_std_byday ;
after_success_std_actual_shuffle_diff_overallmouse = mean(after_success_std_actual_shuffle_diff,2);

%Now we want to see how many times the std in actual is less than that
%in the shuffled
success_std_v_shuffled_pval=[];
after_success_duration_std_inline = reshape(after_success_duration_std,[size(after_success_duration_std,1)*size(after_success_duration_std,2),1]);
for perm_it = 1:length(after_success_duration_std_inline)
    success_std_v_shuffled_pval =[success_std_v_shuffled_pval; sum(std_Shuffled_Length_After_Success_All_Shuffles(perm_it,:) <  after_success_duration_std_inline(perm_it))/size(std_Shuffled_Length_After_Success_All_Shuffles,2)];
end
success_std_v_shuffled_pval_bymouseday = reshape(success_std_v_shuffled_pval,[total_mice,day_number]);
%how many days out of the total number is the SD in actual smaller than the
%shuffled version?
sum(success_std_v_shuffled_pval_bymouseday <= .05,2)/day_number;

%% Mean consec and prop in upstate vs shuffled
Up_State_Consecutive_Mean_All = reshape(Up_State_Consecutive_Mean_All,[total_mice,day_number]);
Up_state_prop_All= reshape(Up_state_prop_All,[total_mice,day_number]);
Up_State_Consecutive_Mean_All_Shuf_dist_mean =reshape(Up_State_Consecutive_Mean_All_Shuf_dist_mean,[total_mice,day_number]);
Up_state_prop_All_Shuf_dist_mean= reshape(Up_state_prop_All_Shuf_dist_mean,[total_mice,day_number]);

%% loop through the cell to plop all animals together in a single matrix, rather than a cell.
%we will use this long matrices to create a table for LMEs
LP_Durations_All = [];
Logical_LP_All =[];
LP_Timestamps_All =[];
Rewarded_Durations_All = [];
Unrewarded_Durations_All =[];
n_minus_one_Rewarded_Durations_All = [];
n_minus_two_Rewarded_Durations_All =[];
n_minus_three_Rewarded_Durations_All =[];
n_minus_four_Rewarded_Durations_All =[];
n_minus_five_Rewarded_Durations_All =[];
n_minus_six_Rewarded_Durations_All =[];
n_minus_seven_Rewarded_Durations_All =[];
n_minus_eight_Rewarded_Durations_All =[];
n_minus_nine_Rewarded_Durations_All =[];
n_minus_ten_Rewarded_Durations_All =[];
IPI_All =[];
Shuffled_Durations_All =[];
Shufffled_Logical_Lever_Press_All =[];
HE_Indicator_All =[];
HE_n_1_Indicator_All = []; 
HE_n_2_Indicator_All = []; 
HE_n_3_Indicator_All = [];
HE_n_4_Indicator_All = [];
HE_n_5_Indicator_All = [];
HE_n_6_Indicator_All = [];
criteria_indicator_All =[];
criteria_percent_indicator_All =[];
moving_mean_All=[];
uppersum_sem_All =[]; 
up_state_idx_All=[]; 
up_state_idx_n1_All=[];
up_state_idx_n2_All=[];
up_state_idx_n3_All=[];
up_state_idx_n4_All=[];
up_state_idx_n5_All=[];
up_state_idx_n6_All=[];
avg_duration_indicator_All =[];
reward_indicator_All =[];
reward_indicator_n1_All =[];
reward_indicator_n2_All =[];
reward_indicator_n3_All =[];
reward_indicator_n4_All =[];
moving_average_lp_length_n7andback_All =[];

stim_indicator_All=[];
stim_indicator_n1_All=[];
stim_indicator_n2_All=[];
stim_indicator_n3_All=[];
stim_indicator_n4_All=[];
stim_indicator_n5_All=[];
stim_indicator_n6_All=[];

histogram_counts_All =[];
histogram_bins_All=[];
histogam_indexs_All=[];       
%individual LME predictions
correctCIprop_criteria_All =[];
correctCIprop_criteria_complex_All =[];
indiv_lme_10_indiv_se_All =[];
indiv_lme_10_indiv_coef_All =[];
indiv_lme_reduced_se_All =[];
indiv_lme_reduced_coef_All =[];
r_squared_criteria_All=[];
r_squared_adjusted_criteria_All = [];
r_squared_complex_criteria_All = [];
r_squared_adjusted_complex_criteria_All = [];

for i = 1:length(LP_Durations_Cell)    
     correctCIprop_criteria_All =[correctCIprop_criteria_All; indiv_correctCIprop_criteria_Cell{i}];
     correctCIprop_criteria_complex_All =[correctCIprop_criteria_complex_All; indiv_correctCIprop_criteria_complex_Cell{i}];
     indiv_lme_10_indiv_se_All =[indiv_lme_10_indiv_se_All indiv_lme_10_indiv_se_Cell{i}];
     indiv_lme_10_indiv_coef_All =[indiv_lme_10_indiv_coef_All indiv_lme_10_indiv_coef_Cell{i}];
     indiv_lme_reduced_se_All =[indiv_lme_reduced_se_All indiv_lme_reduced_se_Cell{i}];
     indiv_lme_reduced_coef_All =[indiv_lme_reduced_coef_All indiv_lme_reduced_coef_Cell{i}];
     r_squared_criteria_All=[r_squared_criteria_All; r_squared_criteria_Cell{i}];
     r_squared_adjusted_criteria_All = [r_squared_adjusted_criteria_All; r_squared_adjusted_criteria_Cell{i}];
     r_squared_complex_criteria_All = [r_squared_complex_criteria_All; r_squared_complex_criteria_Cell{i}];
     r_squared_adjusted_complex_criteria_All = [r_squared_adjusted_complex_criteria_All; r_squared_adjusted_complex_criteria_Cell{i}];
     LP_Durations_All = [LP_Durations_All; LP_Durations_Cell{i}]; 
     Logical_LP_All = [Logical_LP_All; Logcial_LP_Cell{i}];
     LP_Timestamps_All = [LP_Timestamps_All; LP_Timestamps_Cell{i}];
     Rewarded_Durations_All = [Rewarded_Durations_All; Rewarded_Durations_Cell{i}];
     Unrewarded_Durations_All = [Unrewarded_Durations_All; Unrewarded_Durations_Cell{i}];
     IPI_All =[IPI_All; IPI_Cell{i}];
     Shuffled_Durations_All = [Shuffled_Durations_All; Shuffled_LP_Duration_Cell{i}];
     Shufffled_Logical_Lever_Press_All =[Shuffled_Durations_All; Shufffled_Logical_Lever_Press_Cell{i}];
     HE_Indicator_All = [HE_Indicator_All; HE_Indicator_Cell{i}];
     HE_n_1_Indicator_All =[HE_n_1_Indicator_All; HE_n_1_Indicator_Cell{i}]; 
     HE_n_2_Indicator_All = [HE_n_2_Indicator_All; HE_n_2_Indicator_Cell{i}]; 
     HE_n_3_Indicator_All = [HE_n_3_Indicator_All; HE_n_3_Indicator_Cell{i}];
     HE_n_4_Indicator_All = [HE_n_4_Indicator_All; HE_n_4_Indicator_Cell{i}];
     HE_n_5_Indicator_All = [HE_n_5_Indicator_All; HE_n_5_Indicator_Cell{i}];
     HE_n_6_Indicator_All = [HE_n_6_Indicator_All; HE_n_6_Indicator_Cell{i}];
     criteria_indicator_All =[criteria_indicator_All;criteria_indicator_Cell{i}];
     criteria_percent_indicator_All =[criteria_percent_indicator_All;criteria_percent_indicator_Cell{i}];
     avg_duration_indicator_All = [avg_duration_indicator_All; avg_duration_indicator_Cell{i}];
     moving_mean_All = [moving_mean_All; moving_mean_cell{i}];
     up_state_idx_n1_All =[up_state_idx_n1_All; up_state_idx_n1_Cell{i}];
     up_state_idx_n2_All =[up_state_idx_n2_All; up_state_idx_n2_Cell{i}];
     up_state_idx_n3_All =[up_state_idx_n3_All; up_state_idx_n3_Cell{i}];
     up_state_idx_n4_All =[up_state_idx_n4_All; up_state_idx_n4_Cell{i}];
     up_state_idx_n5_All =[up_state_idx_n5_All; up_state_idx_n5_Cell{i}];
     up_state_idx_n6_All =[up_state_idx_n6_All; up_state_idx_n6_Cell{i}];
     up_state_idx_All =[up_state_idx_All; up_state_idx_cell{i}];
     uppersum_sem_All =[uppersum_sem_All; uppersum_sem_cell{i}];
     moving_average_lp_length_n7andback_All =[moving_average_lp_length_n7andback_All; moving_average_lp_length_n7andback_Cell{i}];
     
     stim_indicator_All = [stim_indicator_All; stim_indicator_cell{i}];
     stim_indicator_n1_All = [stim_indicator_n1_All; stim_indicator_n1_cell{i}];
     stim_indicator_n2_All = [stim_indicator_n2_All; stim_indicator_n2_cell{i}];
     stim_indicator_n3_All = [stim_indicator_n3_All; stim_indicator_n3_cell{i}];
     stim_indicator_n4_All = [stim_indicator_n4_All; stim_indicator_n4_cell{i}];
     stim_indicator_n5_All = [stim_indicator_n5_All; stim_indicator_n5_cell{i}];
     stim_indicator_n6_All = [stim_indicator_n6_All; stim_indicator_n6_cell{i}];
         
     reward_indicator_All = [reward_indicator_All; reward_indicator_Cell{i}];
     reward_indicator_n1_All =[reward_indicator_n1_All; reward_indicator_n1_Cell{i}];
     reward_indicator_n2_All =[reward_indicator_n2_All; reward_indicator_n2_Cell{i}];
     reward_indicator_n3_All =[reward_indicator_n3_All; reward_indicator_n3_Cell{i}];
     reward_indicator_n4_All =[reward_indicator_n4_All; reward_indicator_n4_Cell{i}];
     histogram_counts_All =[histogram_counts_All;  histogram_counts_Cell{i}];
     histogram_bins_All=[histogram_bins_All; histogram_bins_Cell{i}];
     histogam_indexs_All=[histogam_indexs_All; histogam_indexs_Cell{i}];      
end
  
%% the n-back variables are strucuted a bit differently, so loop through
%those on their own
press_n_back_Durations_Cell ={};
logical_n_back_cell ={};
ipi_n_back_Durations_Cell ={};
for i = 1:day_number
     ipi_n_back_Durations_Cell = [ipi_n_back_Durations_Cell Data.Day(i).Mouse(1:total_mice).n_back_IPIs];
     press_n_back_Durations_Cell = [press_n_back_Durations_Cell Data.Day(i).Mouse(1:total_mice).n_back_Lengths];
     logical_n_back_cell = [logical_n_back_cell Data.Day(i).Mouse(1:total_mice).Logical_n_back_Lengths];
end

% Create variables for n-back durations, outcomes, ipis, etc.
n_back_All = {};
logical_n_back_All ={};
ipi_n_back_All ={};
    for i = 1:10
    n_back_All = [n_back_All; press_n_back_Durations_Cell{i,1:end}];
    logical_n_back_All =[logical_n_back_All; logical_n_back_cell{i,1:end}];
    ipi_n_back_All =[ipi_n_back_All; ipi_n_back_Durations_Cell{i,1:end}];
    end
   
%end result is a matrix where rows are n back durations, columns are
%lever presses (lined up so n-1 and n-2 are from the same n lever press
%across all animals/days) 
   n_back_All_matrix = cell2mat( n_back_All);
   %check to make sure the cell array is the same length for all nbacks.
   %Need to remove it if animals made fewer than 10 lps, as this will cause
   %the cell to have unequal size and be unable to convert. 
   check_length = length(logical_n_back_All{1,1});
   
   for check_length_it = 1:length(logical_n_back_All)
       if length(logical_n_back_All{check_length_it}) > check_length_it
           logical_n_back_All{check_length_it} =  logical_n_back_All{check_length_it}(1:check_length);
       end
   end           
   
logical_n_back_All_matrix = cell2mat(logical_n_back_All);
ipi_n_back_matrix = cell2mat(ipi_n_back_All);

 %add in press n duration to the n-back matrix
LP_n_back_DataRaster = [n_back_All_matrix; LP_Durations_All';];
%and press n reward logical
logical_LP_n_back_DataRaster = [logical_n_back_All_matrix; Logical_LP_All';];

%regressions require variables in a table, with each column being a
%variable, so flip
LP_n_back_DataRaster = LP_n_back_DataRaster';
logical_LP_n_back_DataRaster = logical_LP_n_back_DataRaster';
ipi_n_back_matrix = ipi_n_back_matrix';

%% Shuffled order Variables
%we need to do the same thing as above, but for shuffled versions of the
%variables
%this is slightly different, since we will have as many versions of the
%shuffled variables as we had shuffle iterations in MEDPC_Behavior_Extract_For_Regressions
%the cells will be number of mouse*days long, and each cell will have an
%individual animal's data, including all shuffle iterations
shuf_n_back_Lengths_distribution_Cell ={};
Shuffled_Durs_Distribution_Cell ={};
shuf_Logical_n_back_Lengths_distribution_Cell ={};
shuffled_ipi1_Cell = {};
shuffled_ipi2_Cell = {};
shuffled_ipi3_Cell = {};
shuffled_he1_Cell ={};
shuffled_he2_Cell ={};
shuffled_he3_Cell ={};
shuffled_he4_Cell ={};
shuffled_MA_Cell ={};
shuffled_ts_Cell ={};
shuffled_upstate_n1_Cell ={};
shuffled_total_reward_indicator_Cell={};
shuffled_MA_n7_Distribution_Cell={};

Shuffled_stim_ind_all_Cell = {};
Shuffled_stim_ind_all_n1_Cell={};
Shuffled_stim_ind_all_n2_Cell={};
Shuffled_stim_ind_all_n3_Cell={};
Shuffled_stim_ind_all_n4_Cell = {};
Shuffled_stim_ind_all_n5_Cell = {};
Shuffled_stim_ind_all_n6_Cell = {};

for i = 1:day_number
    shuf_n_back_Lengths_distribution_Cell = [shuf_n_back_Lengths_distribution_Cell Data.Day(i).Mouse(1:total_mice).shuf_n_back_Lengths_distribution];
    shuf_Logical_n_back_Lengths_distribution_Cell = [shuf_Logical_n_back_Lengths_distribution_Cell Data.Day(i).Mouse(1:total_mice).shuf_Logical_n_back_Lengths_distribution];
    Shuffled_Durs_Distribution_Cell = [Shuffled_Durs_Distribution_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_Durs_Distribtuion];
    shuffled_ipi1_Cell = [shuffled_ipi1_Cell Data.Day(i).Mouse(1:total_mice).shuffled_ipi1_Distribution];
    shuffled_ipi2_Cell = [shuffled_ipi2_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_ipi2_Distribution];
    shuffled_ipi3_Cell = [shuffled_ipi3_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_ipi3_Distribution];
    shuffled_he1_Cell = [shuffled_he1_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he1_Distribution];
    shuffled_he2_Cell = [shuffled_he2_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he2_Distribution];
    shuffled_he3_Cell = [shuffled_he3_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he3_Distribution];
    shuffled_he4_Cell = [shuffled_he4_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_he4_Distribution];
    shuffled_MA_Cell = [shuffled_MA_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_MA_Distribution];
    shuffled_ts_Cell = [shuffled_ts_Cell  Data.Day(i).Mouse(1:total_mice).Shuffled_ts_Distribution];
    shuffled_upstate_n1_Cell =[shuffled_upstate_n1_Cell  Data.Day(i).Mouse(1:total_mice).Shuffled_Up_State_idx_n1_Distribution];
    shuffled_total_reward_indicator_Cell=[shuffled_total_reward_indicator_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_total_reward_indicator_dist];
    shuffled_MA_n7_Distribution_Cell=[shuffled_MA_n7_Distribution_Cell  Data.Day(i).Mouse(1:total_mice).shuffled_MA_n7_Distribution ];
    
    Shuffled_stim_ind_all_Cell = [Shuffled_stim_ind_all_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_Distribution];
    Shuffled_stim_ind_all_n1_Cell = [Shuffled_stim_ind_all_n1_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n1_Distribution];
    Shuffled_stim_ind_all_n2_Cell = [Shuffled_stim_ind_all_n2_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n2_Distribution];
    Shuffled_stim_ind_all_n3_Cell = [Shuffled_stim_ind_all_n3_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n3_Distribution];
    Shuffled_stim_ind_all_n4_Cell = [Shuffled_stim_ind_all_n4_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n4_Distribution];
    Shuffled_stim_ind_all_n5_Cell = [Shuffled_stim_ind_all_n5_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n5_Distribution];
    Shuffled_stim_ind_all_n6_Cell = [Shuffled_stim_ind_all_n6_Cell Data.Day(i).Mouse(1:total_mice).Shuffled_stim_ind_all_n5_Distribution];  
end

%now moving to the "all" versions of the cell, columns will be invidiual
%shuffle iterations, while all the mouse/day data will be lumped together
%in the rows, just as in the non-shuffled variables
shuffled_Durs_Dist_All = [];
shuffled_ipi1_All = [];
shuffled_ipi2_All = [];
shuffled_ipi3_All = [];
shuffled_he1_All =[];
shuffled_he2_All =[];
shuffled_he3_All =[];
shuffled_he4_All =[];
shuffled_MA_All =[];
shuffled_ts_All =[];
shuffled_upstate_n1_All=[];
shuffled_total_reward_indicator_All=[];
shuffled_MA_n7_Distribution_All=[];

    Shuffled_stim_ind_all_All = [];
    Shuffled_stim_ind_all_n1_All = [];
    Shuffled_stim_ind_all_n2_All = [];
    Shuffled_stim_ind_all_n3_All = [];
    Shuffled_stim_ind_all_n4_All = [];
    Shuffled_stim_ind_all_n5_All = [];
    Shuffled_stim_ind_all_n6_All = [];

       
for i = 1:size(Shuffled_Durs_Distribution_Cell,2)
    shuffled_Durs_Dist_All = [shuffled_Durs_Dist_All; Shuffled_Durs_Distribution_Cell{1,i}]; 
    shuffled_ipi1_All = [shuffled_ipi1_All; shuffled_ipi1_Cell{1,i}];
    shuffled_ipi2_All = [shuffled_ipi2_All; shuffled_ipi2_Cell{1,i}];
    shuffled_ipi3_All = [shuffled_ipi3_All; shuffled_ipi3_Cell{1,i}];
    shuffled_he1_All = [shuffled_he1_All; shuffled_he1_Cell{1,i}];
    shuffled_he2_All = [shuffled_he2_All; shuffled_he2_Cell{1,i}];
    shuffled_he3_All = [shuffled_he3_All; shuffled_he3_Cell{1,i}];
    shuffled_he4_All = [shuffled_he4_All; shuffled_he4_Cell{1,i}];
    shuffled_MA_All = [shuffled_MA_All; shuffled_MA_Cell{1,i}];
    shuffled_ts_All = [shuffled_ts_All; shuffled_ts_Cell{1,i}];
    shuffled_upstate_n1_All = [shuffled_upstate_n1_All; shuffled_upstate_n1_Cell{1,i}];
    shuffled_total_reward_indicator_All=[shuffled_total_reward_indicator_All; shuffled_total_reward_indicator_Cell{1,i}];
    shuffled_MA_n7_Distribution_All=[shuffled_MA_n7_Distribution_All;  shuffled_MA_n7_Distribution_Cell{1,i}];
    Shuffled_stim_ind_all_All = [Shuffled_stim_ind_all_All; Shuffled_stim_ind_all_Cell{1,i}];
    Shuffled_stim_ind_all_n1_All = [Shuffled_stim_ind_all_n1_All; Shuffled_stim_ind_all_n1_Cell{1,i}];
    Shuffled_stim_ind_all_n2_All = [Shuffled_stim_ind_all_n2_All; Shuffled_stim_ind_all_n2_Cell{1,i}];
    Shuffled_stim_ind_all_n3_All = [Shuffled_stim_ind_all_n3_All; Shuffled_stim_ind_all_n3_Cell{1,i}];
    Shuffled_stim_ind_all_n4_All = [Shuffled_stim_ind_all_n4_All; Shuffled_stim_ind_all_n4_Cell{1,i}];
    Shuffled_stim_ind_all_n5_All = [Shuffled_stim_ind_all_n5_All; Shuffled_stim_ind_all_n5_Cell{1,i}];
    Shuffled_stim_ind_all_n6_All = [Shuffled_stim_ind_all_n6_All; Shuffled_stim_ind_all_n6_Cell{1,i}];
end

%% Now we need to create a table that includes all our variables,
%each variable to a column, rows are all the data together across mice/days
% logical for rewards
logical_LP_10_back_DataRaster = logical_LP_n_back_DataRaster(:,1:10);
logical_LP_10_back_DataRaster = [logical_LP_10_back_DataRaster logical_LP_n_back_DataRaster(:,11)];
%add in the durations as well
Logical_and_Durations  = [logical_LP_10_back_DataRaster LP_n_back_DataRaster(:,1:10)];
Logical_and_Durations = [Logical_and_Durations LP_n_back_DataRaster(:,11)];
%convert to a table, name all the variables. 

%n_minus_x_All = Logical for reward
%n_minus_one_Durations_All = actual durations in ms
T10_Logical_and_Continuous = array2table(Logical_and_Durations,'VariableNames',{'n_minus_one_All',...
    'n_minus_two_All', 'n_minus_three_All', 'n_minus_four_All',...
    'n_minus_five_All', 'n_minus_six_All', 'n_minus_seven_All',...
    'n_minus_eight_All', 'n_minus_nine_All', 'n_minus_ten_All','LP_All','n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All', 'n_minus_four_Durations_All',...
    'n_minus_five_Durations_All', 'n_minus_six_Durations_All', 'n_minus_seven_Durations_All',...
    'n_minus_eight_Durations_All', 'n_minus_nine_Durations_All', 'n_minus_ten_Durations_All','LP_Durations_All'});

%add in indicator variables
T10_Logical_and_Continuous.mouse_indicator =mouse_indicator; %each lp tagged with mouse number
T10_Logical_and_Continuous.day_indicator = day_indicator; %each lp tagged with session number
%make sure logical variables are coded as categorical
T10_Logical_and_Continuous.mouse_indicator = categorical(T10_Logical_and_Continuous.mouse_indicator);
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
%indicator of mouse's avg. criteria% on a given day
T10_Logical_and_Continuous.criteria_percent_indicator = criteria_percent_indicator_All;
%Did mice make a HE on n-back presses?
T10_Logical_and_Continuous.HE_n_1_Indicator_All = HE_n_1_Indicator_All;
T10_Logical_and_Continuous.HE_n_1_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_1_Indicator_All);
T10_Logical_and_Continuous.HE_n_2_Indicator_All = HE_n_2_Indicator_All;
T10_Logical_and_Continuous.HE_n_2_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_2_Indicator_All);
T10_Logical_and_Continuous.HE_n_3_Indicator_All = HE_n_3_Indicator_All;
T10_Logical_and_Continuous.HE_n_3_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_3_Indicator_All);
T10_Logical_and_Continuous.HE_n_4_Indicator_All = HE_n_4_Indicator_All;
T10_Logical_and_Continuous.HE_n_4_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_4_Indicator_All);
T10_Logical_and_Continuous.HE_n_5_Indicator_All = HE_n_5_Indicator_All;
T10_Logical_and_Continuous.HE_n_5_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_5_Indicator_All);
T10_Logical_and_Continuous.HE_n_6_Indicator_All = HE_n_6_Indicator_All;
T10_Logical_and_Continuous.HE_n_6_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_6_Indicator_All);
%timestamps of when LPs happened 
T10_Logical_and_Continuous.LP_Timestamps_All = LP_Timestamps_All;
%add in IPI between n and various n-backs
T10_Logical_and_Continuous.ipi1 = ipi_n_back_matrix(:,1);
T10_Logical_and_Continuous.ipi2 = ipi_n_back_matrix(:,2);
T10_Logical_and_Continuous.ipi3 = ipi_n_back_matrix(:,3);
T10_Logical_and_Continuous.ipi4 = ipi_n_back_matrix(:,4);
T10_Logical_and_Continuous.ipi5 = ipi_n_back_matrix(:,5);
T10_Logical_and_Continuous.ipi6 = ipi_n_back_matrix(:,6);

%was n-1 in an upstate (from CUSUM)
T10_Logical_and_Continuous.up_state_idx_n1_All =up_state_idx_n1_All;
T10_Logical_and_Continuous.up_state_idx_n1_All = categorical(T10_Logical_and_Continuous.up_state_idx_n1_All);
T10_Logical_and_Continuous.up_state_idx_n2_All =up_state_idx_n2_All;
T10_Logical_and_Continuous.up_state_idx_n2_All = categorical(T10_Logical_and_Continuous.up_state_idx_n2_All);
T10_Logical_and_Continuous.up_state_idx_n3_All =up_state_idx_n3_All;
T10_Logical_and_Continuous.up_state_idx_n3_All = categorical(T10_Logical_and_Continuous.up_state_idx_n3_All);
T10_Logical_and_Continuous.up_state_idx_n4_All =up_state_idx_n4_All;
T10_Logical_and_Continuous.up_state_idx_n4_All = categorical(T10_Logical_and_Continuous.up_state_idx_n4_All);
T10_Logical_and_Continuous.up_state_idx_n5_All =up_state_idx_n5_All;
T10_Logical_and_Continuous.up_state_idx_n5_All = categorical(T10_Logical_and_Continuous.up_state_idx_n5_All);
T10_Logical_and_Continuous.up_state_idx_n6_All =up_state_idx_n6_All;
T10_Logical_and_Continuous.up_state_idx_n6_All = categorical(T10_Logical_and_Continuous.up_state_idx_n6_All);

%avg duration on a day
T10_Logical_and_Continuous.avg_duration_indicator = avg_duration_indicator_All;
%moving avg duration
T10_Logical_and_Continuous.moving_average_lp_length_n7andback_All=moving_average_lp_length_n7andback_All;
%LP_All is logical for reward, and then n-back logicals
T10_Logical_and_Continuous.LP_All =categorical(T10_Logical_and_Continuous.LP_All);
T10_Logical_and_Continuous.n_minus_one_All =categorical(T10_Logical_and_Continuous.n_minus_one_All);
T10_Logical_and_Continuous.n_minus_two_All =categorical(T10_Logical_and_Continuous.n_minus_two_All);
T10_Logical_and_Continuous.n_minus_three_All =categorical(T10_Logical_and_Continuous.n_minus_three_All);
T10_Logical_and_Continuous.n_minus_four_All =categorical(T10_Logical_and_Continuous.n_minus_four_All);
T10_Logical_and_Continuous.n_minus_five_All =categorical(T10_Logical_and_Continuous.n_minus_five_All);
T10_Logical_and_Continuous.n_minus_six_All =categorical(T10_Logical_and_Continuous.n_minus_six_All);

%indicator for what the hold down criteria was (e.g., 800 or 1600)
T10_Logical_and_Continuous.criteria_indicator_All = criteria_indicator_All;
T10_Logical_and_Continuous.criteria_indicator_All = categorical(T10_Logical_and_Continuous.criteria_indicator_All);
% logical indicators for if stimulation occurred
T10_Logical_and_Continuous.stim_indicator_All = stim_indicator_All;
T10_Logical_and_Continuous.stim_indicator_n1_All = stim_indicator_n1_All;
T10_Logical_and_Continuous.stim_indicator_n2_All = stim_indicator_n2_All;
T10_Logical_and_Continuous.stim_indicator_n3_All = stim_indicator_n3_All;
T10_Logical_and_Continuous.stim_indicator_n4_All = stim_indicator_n4_All;
T10_Logical_and_Continuous.stim_indicator_n5_All = stim_indicator_n5_All;
T10_Logical_and_Continuous.stim_indicator_n6_All = stim_indicator_n6_All;

%logical indicators of Reward - these only differ from the n_minus_x_All
%variables when reward is probabalistic.
T10_Logical_and_Continuous.reward_indicator_All = reward_indicator_All;
T10_Logical_and_Continuous.reward_indicator_n1_All =reward_indicator_n1_All;
T10_Logical_and_Continuous.reward_indicator_n2_All =reward_indicator_n2_All; 
T10_Logical_and_Continuous.reward_indicator_n3_All =reward_indicator_n3_All;
T10_Logical_and_Continuous.reward_indicator_n4_All =reward_indicator_n4_All; 
T10_Logical_and_Continuous.reward_indicator_n1_All = categorical(T10_Logical_and_Continuous.reward_indicator_n1_All);
T10_Logical_and_Continuous.reward_indicator_n2_All = categorical(T10_Logical_and_Continuous.reward_indicator_n2_All);
T10_Logical_and_Continuous.reward_indicator_n3_All = categorical(T10_Logical_and_Continuous.reward_indicator_n3_All);
T10_Logical_and_Continuous.reward_indicator_n4_All = categorical(T10_Logical_and_Continuous.reward_indicator_n4_All);

%% optional remove Nans for export to R
% T10_nansgone = T10_Logical_and_Continuous;
% T10_nansgone =T10_nansgone(~any(ismissing(T10_nansgone),2),:);
% writetable(T10_nansgone,'33-3-all-time.csv','QuoteStrings',true) %save table as
% csv
 %% This section is primarily applicable for optogenetic Stimulation, and is used to find variables relative to stimulation
 %e.g., press duration during/after stimulation to look for direct effect
 %also some calculations of overall mouse mean durations and such
 %get lps and nbacks aligned with stim logical (1 = stim, 0 = no stim)
 no_stim_idx =  find(T10_Logical_and_Continuous.stim_indicator_All(:) ~=1);
 stim_idx =find(T10_Logical_and_Continuous.stim_indicator_All(:)==1);
 %press index after stim/nostim
 no_stim_idx_n1 =  find(T10_Logical_and_Continuous.stim_indicator_n1_All(:) ~=1);
 stim_idx_n1 =find(T10_Logical_and_Continuous.stim_indicator_n1_All(:)==1); 
%press durations 
 lp_stim = T10_Logical_and_Continuous.LP_Durations_All(stim_idx);
 lp_no_stim = T10_Logical_and_Continuous.LP_Durations_All(no_stim_idx);
 %n + 1 press durtions
 lp_stim_n1 = T10_Logical_and_Continuous.LP_Durations_All(stim_idx_n1);
 lp_no_stim_n1 = T10_Logical_and_Continuous.LP_Durations_All(no_stim_idx_n1);

 %ipi after stim
 ipi_stim_n1 = T10_Logical_and_Continuous.ipi1(stim_idx_n1);
 ipi_no_stim_n1 = T10_Logical_and_Continuous.ipi1(no_stim_idx_n1);
 
 %Get Proportion of HEs after stim or after nostim
 he_stim_n1 = T10_Logical_and_Continuous.HE_n_1_Indicator_All(stim_idx_n1);
 he_stim_n1=double(he_stim_n1);
 sum(he_stim_n1(he_stim_n1==2))/length(he_stim_n1);
 he_no_stim_n1 = T10_Logical_and_Continuous.HE_n_1_Indicator_All(no_stim_idx_n1);
 he_no_stim_n1=double(he_no_stim_n1);
 sum(he_no_stim_n1(he_no_stim_n1==2))/length(he_no_stim_n1);

%get means per animal
%briefly convert categorical mouse/day indicators to double to use indexing
T10_Logical_and_Continuous.day_indicator = double(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =double(T10_Logical_and_Continuous.mouse_indicator);
%initialize variables
mouse_means =[];
mouse_std =[];
mouse_mean_over_std =[];
mouse_med=[];
mouse_iqr=[];
mouse_iqr_over_med=[];
mouse_no_stim_mean =[];
mouse_stim_mean =[];
mouse_stim_n1_mean =[];
mouse_no_stim_n1_mean =[];
mouse_stim_n1_ipi_mean =[];
mouse_no_stim_n1_ipi_mean =[];
mouse_stim_n1_std =[];
mouse_no_stim_n1_std =[];
mouse_stim_n_and_n1 =[];
mouse_no_stim_n_and_n1=[];

%create day and mouse iteration variables to loop through
it_d =length(unique(T10_Logical_and_Continuous.day_indicator));
it_m = length(unique(T10_Logical_and_Continuous.mouse_indicator));

for day_mean = 1:it_d 
    %loop through days uisng da_index
     day_idx = find(T10_Logical_and_Continuous.day_indicator ==day_mean);
   
for mouse_mean = 1:it_m
    %then within each day, loop through mice, using the index
    mouse_idx = find(T10_Logical_and_Continuous.mouse_indicator ==mouse_mean);
    %the intersection of mouse/day will give you data for individual mouse
    %sessions
    mouse_by_day_idx = intersect(day_idx,mouse_idx);
    
    %now calculate mean, std, med, interquartile range (IQR) overall per
    %mouse/day
    mouse_means = [mouse_means mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    mouse_std = [mouse_std std(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    mouse_mean_over_std =[mouse_mean_over_std mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))/std(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    mouse_med=[mouse_med median(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    mouse_iqr=[mouse_iqr iqr(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    mouse_iqr_over_med=[mouse_iqr_over_med iqr(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))/ median(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx))];
    
    %and now use the intersection of stim/nostim index to find lps that
    %were/not stimulated
    mouse_by_day_stim_idx = intersect(mouse_by_day_idx,stim_idx);
    mouse_by_day_no_stim_idx = intersect(mouse_by_day_idx,no_stim_idx);
    %get means, ipi, etc. for press with stimulation/w/o
    mouse_stim_mean = [mouse_stim_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx))];
    mouse_no_stim_mean = [mouse_no_stim_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx))];
    
    %now for n+1
    mouse_by_day_stim_idx_n1 = intersect(mouse_by_day_idx,stim_idx_n1);
    mouse_by_day_no_stim_idx_n1 = intersect(mouse_by_day_idx,no_stim_idx_n1);
    
    mouse_stim_n1_mean = [mouse_stim_n1_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1))];
    mouse_no_stim_n1_mean = [mouse_no_stim_n1_mean mean(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1))];
    mouse_stim_n1_ipi_mean =[mouse_stim_n1_ipi_mean nanmean(T10_Logical_and_Continuous.ipi1(mouse_by_day_stim_idx_n1))];
    mouse_no_stim_n1_ipi_mean =[mouse_no_stim_n1_ipi_mean nanmean(T10_Logical_and_Continuous.ipi1(mouse_by_day_no_stim_idx_n1))];
    mouse_stim_n1_std = [mouse_stim_n1_std nanstd(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1))];
    mouse_no_stim_n1_std = [mouse_no_stim_n1_std nanstd(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1))];
    
%we want to match the number of no stim to stim trials to calculate
%correlations. E.g., how related are press n and n - 1 normally - does this
%change with stimulation?

%paste n and n+1 together into an array - will have to pad the n1
%with nans appropriately to match length
mouse_stim_lps =T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx);
mouse_stim_n1_lps = T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1);
%use this to match length of n and n+1
if length(mouse_stim_lps) > length(mouse_stim_n1_lps)
    mouse_stim_n1_lps = [mouse_stim_n1_lps;NaN];
end
  
stim_n_and_n1 = [mouse_stim_lps mouse_stim_n1_lps];
%all mice lumped together
mouse_stim_n_and_n1 = [mouse_stim_n_and_n1; stim_n_and_n1];

%repeat for no stim.
mouse_no_stim_lps =T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx);
mouse_no_stim_n1_lps = T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1);
mouse_no_stim_n1_lps(1) =[];
mouse_no_stim_n1_lps = [mouse_no_stim_n1_lps;NaN];

if length(mouse_no_stim_lps) > length(mouse_no_stim_n1_lps)
    mouse_no_stim_n1_lps = [mouse_no_stim_n1_lps;NaN];
end
  
if length(mouse_no_stim_lps) < length(mouse_no_stim_n1_lps)
    mouse_no_stim_lps = [mouse_no_stim_lps;NaN];
end

no_stim_n_and_n1 = [mouse_no_stim_lps mouse_no_stim_n1_lps];
mouse_no_stim_n_and_n1 = [mouse_no_stim_n_and_n1; no_stim_n_and_n1];
end
end

%reshape to get mouse(row) x day(col), mostly for graphing
mouse_stim_mean = reshape(mouse_stim_mean,[total_mice,day_number]);
mouse_no_stim_mean = reshape(mouse_no_stim_mean,[total_mice,day_number]);
mouse_stim_n1_mean = reshape(mouse_stim_n1_mean,[total_mice,day_number]);
mouse_no_stim_n1_mean = reshape(mouse_no_stim_n1_mean,[total_mice,day_number]);
mouse_stim_n1_ipi_mean = reshape(mouse_stim_n1_ipi_mean,[total_mice,day_number]);
mouse_no_stim_n1_ipi_mean = reshape(mouse_no_stim_n1_ipi_mean,[total_mice,day_number]);
mouse_stim_n1_std = reshape(mouse_stim_n1_std,[total_mice,day_number]);
mouse_no_stim_n1_std = reshape(mouse_no_stim_n1_std,[total_mice,day_number]);
mouse_means = reshape(mouse_means,[total_mice,day_number]);
mouse_std = reshape(mouse_std,[total_mice,day_number]);
mouse_mean_over_std =reshape(mouse_mean_over_std,[total_mice,day_number]);
mouse_med=reshape(mouse_med,[total_mice,day_number]);
mouse_iqr=reshape(mouse_iqr,[total_mice,day_number]);
mouse_iqr_over_med=reshape(mouse_iqr_over_med,[total_mice,day_number]);

%% is there a correlation between individual model fit and task performance?
%export to R and use RMcorr on the table 
table_of_r2 =  array2table(r_squared_criteria_All,'VariableNames',{'R2','Criteria'});
table_of_r2_complex =  array2table(r_squared_complex_criteria_All,'VariableNames',{'R2','Criteria'});
%is there a correlation between n-1 coef value and r2?
table_of_r2.Coef = indiv_lme_10_indiv_coef_All(2,:)';
%or of model predictions
table_of_r2.Predictions = correctCIprop_criteria_All(:,1);
%add in mouse number
table_of_r2.Mouse = mouse_number_for_corrs;
%save to a csv for export to R for use in RMCORR
% writetable(table_of_r2,'m2_dms_sham_diffpreds_simple.csv','QuoteStrings',true) %save table as
%using complex LMEs
table_of_r2_complex.Pred = correctCIprop_criteria_complex_All(:,1);
%add in mouse number
table_of_r2_complex.Mouse = mouse_number_for_corrs;
table_of_r2_complex.Day = day_number_for_corrs;

%save to a csv for export to R for use in RMCORR
% writetable(table_of_r2_complex,'m2_dms_sham_diffpreds_complex.csv','QuoteStrings',true) %save table as
 
%% return the indicator variables to categorical for use in LMEs
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =categorical(T10_Logical_and_Continuous.mouse_indicator);
T10_Logical_and_Continuous.stim_indicator_All = categorical(T10_Logical_and_Continuous.stim_indicator_All);
T10_Logical_and_Continuous.stim_indicator_n1_All = categorical(T10_Logical_and_Continuous.stim_indicator_n1_All);
T10_Logical_and_Continuous.stim_indicator_n2_All = categorical(T10_Logical_and_Continuous.stim_indicator_n2_All);
T10_Logical_and_Continuous.stim_indicator_n3_All = categorical(T10_Logical_and_Continuous.stim_indicator_n3_All);
T10_Logical_and_Continuous.stim_indicator_n4_All = categorical(T10_Logical_and_Continuous.stim_indicator_n4_All);
T10_Logical_and_Continuous.stim_indicator_n5_All = categorical(T10_Logical_and_Continuous.stim_indicator_n5_All);
T10_Logical_and_Continuous.stim_indicator_n6_All = categorical(T10_Logical_and_Continuous.stim_indicator_n6_All);

%% Complex LME Model
%Use BIC selected variables (selected from a full model that included all
%terms and their interactions with durations up to n -6)
model_spec_reduced = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

%create LME, and save the model
lme_reduced = fitlme(T10_Logical_and_Continuous,model_spec_reduced,'DummyVarCoding','reference');
lme_reduced_se = lme_reduced.Coefficients.SE;
lme_reduced_coef =  lme_reduced.Coefficients.Estimate;
lme_reduced_name =  lme_reduced.Coefficients.Name;
lme_reduced_pval =  lme_reduced.Coefficients.pValue;
lme_reduced_t =  lme_reduced.Coefficients.tStat;
lme_reduced_df =  lme_reduced.Coefficients.DF;
lme_reduced_upper =  lme_reduced.Coefficients.Upper;
lme_reduced_lower =  lme_reduced.Coefficients.Lower;
lme_reduced_AIC = lme_reduced.ModelCriterion.AIC;
lme_reduced_anova = anova(lme_reduced);
lme_reduced_fstat = lme_reduced_anova.FStat;
lme_reduced_fpval = lme_reduced_anova.pValue;
[b bnames bstats] =randomEffects(lme_reduced);

%create the same for an lm for easy graphing
model_spec_reduce_lm = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator';
lm_reduced = fitlm(T10_Logical_and_Continuous,model_spec_reduce_lm);

%optional save simple model to fit to other datasets
% save('33_3lm.mat' ,'lm_reduced')

%create n10 back durations model continaing control vars for shuffling
model_spec_10_ma_totrew = 'LP_Durations_All~criteria_percent_indicator+LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +(1|mouse_indicator)+(1|day_indicator)';
lme_10_ma_totrew = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew);
lme_10_ma_totrew_se = lme_10_ma_totrew.Coefficients.SE;
lme_10_ma_totrew_coef =  lme_10_ma_totrew.Coefficients.Estimate;
lme_10_ma_totrew_name =  lme_10_ma_totrew.Coefficients.Name;
lme_10_ma_totrew_pval =  lme_10_ma_totrew.Coefficients.pValue;
lme_10_ma_totrew_t =  lme_10_ma_totrew.Coefficients.tStat;
lme_10_ma_totrew_df =  lme_10_ma_totrew.Coefficients.DF;
lme_10_ma_totrew_upper =  lme_10_ma_totrew.Coefficients.Upper;
lme_10_ma_totrew_lower =  lme_10_ma_totrew.Coefficients.Lower;
lme_10_ma_totrew_AIC = lme_10_ma_totrew.ModelCriterion.AIC;
lme_10_ma_totrew_anova = anova(lme_10_ma_totrew);
lme_10_ma_totrew_fstat =lme_10_ma_totrew_anova.FStat;
lme_10_ma_totrew_fpval = lme_10_ma_totrew_anova.pValue;

%% shuffling
%For the simple LME, we will shuffle each term one at a time (e.g., only
%n-3) and compare the beta coef from the order shuffled datasets to the
%actual coeffs. 
%Also have for the complex model, but only to double check BIC results,
%comment out for speed of running

%initialize variables
%simple n1:n10 and control vars
%the size will be the number of coefficients in model by number of shuffles
shuf_n1_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n1_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n2_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n2_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n3_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n3_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n4_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n4_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n5_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n5_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n6_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n6_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n7_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n7_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n8_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n8_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n9_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n9_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n10_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_n10_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_ts_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_ts_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_crit_simple_coef_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_crit_simple_SE_all = nan(size(lme_10_ma_totrew_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));

%variables from full beh model
shuf_full_ma_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n1dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n2dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n3dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n4dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n5dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n6dur_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_me_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_n1_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_avg_coef_all = nan(size(lme_reduced_coef,1),size(shuf_n_back_Lengths_distribution_Cell,1));

%using AIC instead of Coefs.
shuf_full_n1dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n2dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n3dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n4dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n5dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_n6dur_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ma_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_rew_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ts_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_crit_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_me_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_n1_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_he_avg_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));
shuf_full_ipi2_n2_AIC_all = nan(1,size(shuf_n_back_Lengths_distribution_Cell,1));

%souter loop is for the number of shuffles made 
for j =  1:size(shuf_n_back_Lengths_distribution_Cell,1)
%need to erase the i_shuffled array every time, this will be used to create
%the table for regression
i_shuffled =[];
%inner loop is for the number of mice x sessions to get n-backs
for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
%first shuffle for all mice
i_shuffled = [i_shuffled; shuf_n_back_Lengths_distribution_Cell{j,i}'];
end
%get logical reward n-backs inner loop
shuf_logical_n1=[];
for i = 1:size(shuf_n_back_Lengths_distribution_Cell,2)
    shuf_logical_n1 = [shuf_logical_n1; shuf_Logical_n_back_Lengths_distribution_Cell{j,i}(1,:)'];
end
%both nback durations and logicals for reward
i_shuffled = [i_shuffled shuf_logical_n1];

%add in the appropriate shuffled n lp (ie, the shuffled dist that built the
%n-back array)
  i_shuffled = [i_shuffled shuffled_Durs_Dist_All(:,j)];
 %add in appropriate shuffled moving average and other vars
 i_shuffled = [i_shuffled shuffled_total_reward_indicator_All(:,j)];
 i_shuffled =[i_shuffled shuffled_MA_n7_Distribution_All(:,j)];
 i_shuffled = [i_shuffled shuffled_MA_All(:,j)];
 i_shuffled = [i_shuffled shuffled_ipi1_All(:,j)];
 i_shuffled = [i_shuffled shuffled_ipi2_All(:,j)];
 i_shuffled = [i_shuffled shuffled_ipi3_All(:,j)];
 i_shuffled = [i_shuffled shuffled_he1_All(:,j)];
 i_shuffled = [i_shuffled shuffled_he2_All(:,j)];
 i_shuffled = [i_shuffled shuffled_he3_All(:,j)];
 i_shuffled = [i_shuffled shuffled_he4_All(:,j)];
 i_shuffled = [i_shuffled shuffled_upstate_n1_All(:,j)];
 i_shuffled = [i_shuffled shuffled_ts_All(:,j)];
 
%add in the shuffled indicator variables. As the order of these within day
%cannot be shuffled, they are shuffled across mice/days
shuf_criteria_percent_indicator_All = criteria_percent_indicator_All(randperm(length(criteria_percent_indicator_All)));
shuf_criteria_indicator_All = criteria_indicator_All(randperm(length(criteria_indicator_All)));
shuf_mouse_indicator = mouse_indicator(randperm(length(mouse_indicator)));
shuf_day_indicator = day_indicator(randperm(length(day_indicator)));

  i_shuffled = [i_shuffled shuf_mouse_indicator];
  i_shuffled = [i_shuffled shuf_day_indicator];
  i_shuffled = [i_shuffled shuf_criteria_percent_indicator_All];
  i_shuffled = [i_shuffled shuf_criteria_indicator_All];
  %create a table from the shuffled array
   T_shuffled = array2table(i_shuffled,'VariableNames',{'shuf_n_minus_one_Durations_All',...
    'shuf_n_minus_two_Durations_All', 'shuf_n_minus_three_Durations_All', 'shuf_n_minus_four_Durations_All',...
    'shuf_n_minus_five_Durations_All', 'shuf_n_minus_six_Durations_All', 'shuf_n_minus_seven_Durations_All',...
    'shuf_n_minus_eight_Durations_All', 'shuf_n_minus_nine_Durations_All', 'shuf_n_minus_ten_Durations_All',...
   'shuf_n_minus_one_All','shuf_LP_Durations_All','shuf_total_reward_indicator_All','shuf_moving_average_lp_length_n7andback_All',...
   'shuf_moving_mean','shuf_ipi1','shuf_ipi2','shuf_ipi3','shuf_HE_n_1_Indicator_All',...
   'shuf_HE_n_2_Indicator_All','shuf_HE_n_3_Indicator_All','shuf_HE_n_4_Indicator_All','shuf_up_state_idx_n1_All','shuf_LP_Timestamps_All',...
   'shuf_mouse_indicator','shuf_day_indicator','shuf_criteria_percent_indicator','shuf_criteria_indicator_All'});

%create atable with both the actual and shuffled versions of all variables
T_Both = [T_shuffled T10_Logical_and_Continuous];

%add in stim variables if needed
T_Both.shuf_stim_indicator_All = Shuffled_stim_ind_all_All(:,j);
T_Both.shuf_stim_indicator_n1_All = Shuffled_stim_ind_all_n1_All(:,j);
T_Both.shuf_stim_indicator_n2_All = Shuffled_stim_ind_all_n2_All(:,j);
T_Both.shuf_stim_indicator_n3_All = Shuffled_stim_ind_all_n3_All(:,j);
T_Both.shuf_stim_indicator_n4_All = Shuffled_stim_ind_all_n4_All(:,j);
T_Both.shuf_stim_indicator_n5_All = Shuffled_stim_ind_all_n5_All(:,j);
T_Both.shuf_stim_indicator_n6_All = Shuffled_stim_ind_all_n6_All(:,j);

% shuffle the simple n10 ma tot rew model one at a time for permutation
%tests 
%n1
model_spec_shuf_n1_simple= 'LP_Durations_All~LP_Timestamps_All+shuf_n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n1_simple = fitlme(T_Both,model_spec_shuf_n1_simple);
shuf_n1_simple_coef_all(:,j) =  shuf_n1_simple.Coefficients.Estimate;
shuf_n1_simple_SE_all(:,j) = shuf_n1_simple.Coefficients.SE;
%n2
model_spec_shuf_n2_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + shuf_n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n2_simple = fitlme(T_Both,model_spec_shuf_n2_simple);
shuf_n2_simple_coef_all(:,j) =  shuf_n2_simple.Coefficients.Estimate;
shuf_n2_simple_SE_all(:,j) = shuf_n2_simple.Coefficients.SE;
%n3
model_spec_shuf_n3_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + shuf_n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n3_simple = fitlme(T_Both,model_spec_shuf_n3_simple);
shuf_n3_simple_coef_all(:,j) = shuf_n3_simple.Coefficients.Estimate;;
shuf_n3_simple_SE_all(:,j) = shuf_n3_simple.Coefficients.SE;;
%n4
model_spec_shuf_n4_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + shuf_n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n4_simple = fitlme(T_Both,model_spec_shuf_n4_simple);
shuf_n4_simple_coef_all(:,j) = shuf_n4_simple.Coefficients.Estimate;
shuf_n4_simple_SE_all(:,j) = shuf_n4_simple.Coefficients.SE;
%n5
model_spec_shuf_n5_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+shuf_n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n5_simple = fitlme(T_Both,model_spec_shuf_n5_simple);
shuf_n5_simple_coef_all(:,j) = shuf_n5_simple.Coefficients.Estimate;
shuf_n5_simple_SE_all(:,j) = shuf_n5_simple.Coefficients.SE;
%n6
model_spec_shuf_n6_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+shuf_n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n6_simple = fitlme(T_Both,model_spec_shuf_n6_simple);
shuf_n6_simple_coef_all(:,j) = shuf_n6_simple.Coefficients.Estimate;
shuf_n6_simple_SE_all(:,j) = shuf_n6_simple.Coefficients.SE;
%n7
model_spec_shuf_n7_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+shuf_n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n7_simple = fitlme(T_Both,model_spec_shuf_n7_simple);
shuf_n7_simple_coef_all(:,j) =shuf_n7_simple.Coefficients.Estimate;
shuf_n7_simple_SE_all(:,j) = shuf_n7_simple.Coefficients.SE;
%n8
model_spec_shuf_n8_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+shuf_n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n8_simple = fitlme(T_Both,model_spec_shuf_n8_simple);
shuf_n8_simple_coef_all(:,j) = shuf_n8_simple.Coefficients.Estimate;
shuf_n8_simple_SE_all(:,j) = shuf_n8_simple.Coefficients.SE;
%n9
model_spec_shuf_n9_simple= 'LP_Durations_All~LP_Timestamps_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+shuf_n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n9_simple = fitlme(T_Both,model_spec_shuf_n9_simple);
shuf_n9_simple_coef_all(:,j) = shuf_n9_simple.Coefficients.Estimate;
shuf_n9_simple_SE_all(:,j) = shuf_n9_simple.Coefficients.SE;
%n10
model_spec_shuf_n10_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+shuf_n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_n10_simple = fitlme(T_Both,model_spec_shuf_n10_simple);
shuf_n10_simple_coef_all(:,j) = shuf_n10_simple.Coefficients.Estimate;
shuf_n10_simple_SE_all(:,j) = shuf_n10_simple.Coefficients.SE;
%Timestamps
model_spec_shuf_ts_simple= 'LP_Durations_All~shuf_LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
shuf_ts_simple = fitlme(T_Both,model_spec_shuf_ts_simple);
shuf_ts_simple_coef_all(:,j) = shuf_ts_simple.Coefficients.Estimate;
shuf_ts_simple_SE_all(:,j) = shuf_ts_simple.Coefficients.SE;
%Criteria%
model_spec_shuf_crit_simple= 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +  shuf_criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator)';
shuf_crit_simple = fitlme(T_Both,model_spec_shuf_crit_simple);
shuf_crit_simple_coef_all(:,j) = shuf_crit_simple.Coefficients.Estimate;
shuf_crit_simple_SE_all(:,j) = shuf_crit_simple.Coefficients.SE;

%% shuffle the variables that remained in the  behavioral model selected via BIC for 33-3
% %this is the reduced only, significant stuff model
% % model_spec_reduced = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% % 
% %shuffle variables one at a time 
% %Moving Average Main Effect
% shuf_spec_full_ma = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + shuf_moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ma = fitlme(T_Both,shuf_spec_full_ma);
% shuf_full_ma_coef_all(:,j) =  shuf_full_ma.Coefficients.Estimate;
% shuf_full_ma_AIC_all(:,j) =  shuf_full_ma.ModelCriterion.AIC;
% 
% %now shuffle the n-back durations
% %n1
% shuf_spec_full_n1dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + shuf_n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n1dur = fitlme(T_Both,shuf_spec_full_n1dur);
% shuf_full_n1dur_coef_all(:,j) = shuf_full_n1dur.Coefficients.Estimate;
% shuf_full_n1dur_AIC_all(:,j) =  shuf_full_n1dur.ModelCriterion.AIC;
% 
% shuf_spec_full_n2dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + shuf_n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n2dur = fitlme(T_Both,shuf_spec_full_n2dur);
% shuf_full_n2dur_coef_all(:,j) = shuf_full_n2dur.Coefficients.Estimate;
% shuf_full_n2dur_AIC_all(:,j) =  shuf_full_n2dur.ModelCriterion.AIC;
% 
% shuf_spec_full_n3dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + shuf_n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n3dur = fitlme(T_Both,shuf_spec_full_n3dur);
% shuf_full_n3dur_coef_all(:,j) =  shuf_full_n3dur.Coefficients.Estimate;
% shuf_full_n3dur_AIC_all(:,j) = shuf_full_n3dur.ModelCriterion.AIC;
% 
% shuf_spec_full_n4dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + shuf_n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n4dur = fitlme(T_Both,shuf_spec_full_n4dur);
% shuf_full_n4dur_coef_all(:,j) = shuf_full_n4dur.Coefficients.Estimate;
% shuf_full_n4dur_AIC_all(:,j) = shuf_full_n4dur.ModelCriterion.AIC;
% 
% shuf_spec_full_n5dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + shuf_n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n5dur = fitlme(T_Both,shuf_spec_full_n5dur);
% shuf_full_n5dur_coef_all(:,j) = shuf_full_n5dur.Coefficients.Estimate;
% shuf_full_n5dur_AIC_all(:,j) =  shuf_full_n5dur.ModelCriterion.AIC;
% 
% shuf_spec_full_n6dur = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + shuf_n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_n6dur = fitlme(T_Both,shuf_spec_full_n6dur);
% shuf_full_n6dur_coef_all(:,j) = shuf_full_n6dur.Coefficients.Estimate;
% shuf_full_n6dur_AIC_all(:,j) =  shuf_full_n6dur.ModelCriterion.AIC;
% 
% %% shuffle the interaction variables one at a time
% %ipi me, n1, mov avg
% 
% %ipi1 me
% shuf_spec_full_ipi_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + shuf_ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ipi_me = fitlme(T_Both,shuf_spec_full_ipi_me);
% shuf_full_ipi_me_coef_all(:,j) = shuf_full_ipi_me.Coefficients.Estimate;
% shuf_full_ipi_me_AIC_all(:,j) =  shuf_full_ipi_me.ModelCriterion.AIC;
% %ipi1:n1
% shuf_spec_full_ipi_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  shuf_ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ipi_n1 = fitlme(T_Both,shuf_spec_full_ipi_n1);
% shuf_full_ipi_n1_coef_all(:,j) = shuf_full_ipi_n1.Coefficients.Estimate;
% shuf_full_ipi_n1_AIC_all(:,j) =  shuf_full_ipi_n1.ModelCriterion.AIC;
% %ipi1:MA
% shuf_spec_full_ipi_avg = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  shuf_ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ipi_avg = fitlme(T_Both,shuf_spec_full_ipi_avg);
% shuf_full_ipi_avg_coef_all(:,j) =  shuf_full_ipi_avg.Coefficients.Estimate;
% shuf_full_ipi_avg_AIC_all(:,j) =  shuf_full_ipi_avg.ModelCriterion.AIC;
% 
% %ipi2 me
% shuf_spec_full_ipi2_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   shuf_ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ipi2_me = fitlme(T_Both,shuf_spec_full_ipi2_me);
% shuf_full_ipi2_me_coef_all(:,j) =  shuf_full_ipi2_me.Coefficients.Estimate;
% shuf_full_ipi2_me_AIC_all(:,j) =  shuf_full_ipi2_me.ModelCriterion.AIC;
% %ipi2:n2e
% shuf_spec_full_ipi2_n2 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  shuf_ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ipi2_n2  = fitlme(T_Both,shuf_spec_full_ipi2_n2 );
% shuf_full_ipi2_n2_coef_all(:,j) =  shuf_full_ipi2_n2.Coefficients.Estimate;
% shuf_full_ipi2_n2_AIC_all(:,j) =  shuf_full_ipi2_n2.ModelCriterion.AIC;
% 
% %rew n1 ME
% shuf_spec_full_rew_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + shuf_n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_rew_me = fitlme(T_Both,shuf_spec_full_rew_me);
% shuf_full_rew_me_coef_all(:,j) =  shuf_full_rew_me.Coefficients.Estimate;
% shuf_full_rew_me_AIC_all(:,j) =  shuf_full_rew_me.ModelCriterion.AIC;
% %rew n1:n1
% shuf_spec_full_rew_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:shuf_n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_rew_n1 = fitlme(T_Both,shuf_spec_full_rew_n1);
% shuf_full_rew_n1_coef_all(:,j) =  shuf_full_rew_n1.Coefficients.Estimate;
% shuf_full_rew_n1_AIC_all(:,j) =  shuf_full_rew_n1.ModelCriterion.AIC;
% %rew n1:MA
% shuf_spec_full_rew_avg = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:shuf_n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_rew_avg = fitlme(T_Both,shuf_spec_full_rew_avg);
% shuf_full_rew_avg_coef_all(:,j) =  shuf_full_rew_avg.Coefficients.Estimate;
% shuf_full_rew_avg_AIC_all(:,j) =  shuf_full_rew_avg.ModelCriterion.AIC;
% 
% %TS ME
% shuf_spec_full_ts_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   shuf_LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ts_me = fitlme(T_Both,shuf_spec_full_ts_me);
% shuf_full_ts_me_coef_all(:,j) = shuf_full_ts_me.Coefficients.Estimate;
% shuf_full_ts_me_AIC_all(:,j) =  shuf_full_ts_me.ModelCriterion.AIC;
% %TS:N1
% shuf_spec_full_ts_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:shuf_LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ts_n1 = fitlme(T_Both,shuf_spec_full_ts_n1);
% shuf_full_ts_n1_coef_all(:,j) =  shuf_full_ts_n1.Coefficients.Estimate;
% shuf_full_ts_n1_AIC_all(:,j) =  shuf_full_ts_n1.ModelCriterion.AIC;
% %TS:MA
% shuf_spec_full_ts_avg = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:shuf_LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_ts_avg = fitlme(T_Both,shuf_spec_full_ts_avg);
% shuf_full_ts_avg_coef_all(:,j) =  shuf_full_ts_avg.Coefficients.Estimate;
% shuf_full_ts_avg_AIC_all(:,j) =  shuf_full_ts_avg.ModelCriterion.AIC;
% 
% %Crit% ME
% shuf_spec_full_crit_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + shuf_criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_crit_me = fitlme(T_Both,shuf_spec_full_crit_me);
% shuf_full_crit_me_coef_all(:,j) =  shuf_full_crit_me.Coefficients.Estimate;
% shuf_full_crit_me_AIC_all(:,j) =  shuf_full_crit_me.ModelCriterion.AIC;
% %Crit%:N1
% shuf_spec_full_crit_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:shuf_criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_crit_n1 = fitlme(T_Both,shuf_spec_full_crit_n1);
% shuf_full_crit_n1_coef_all(:,j) =  shuf_full_crit_n1.Coefficients.Estimate;
% shuf_full_crit_n1_AIC_all(:,j) =  shuf_full_crit_n1.ModelCriterion.AIC;
% %Crit%:MA
% shuf_spec_full_crit_avg = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:shuf_criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_crit_avg = fitlme(T_Both,shuf_spec_full_crit_avg);
% shuf_full_crit_avg_coef_all(:,j) =  shuf_full_crit_avg.Coefficients.Estimate;
% shuf_full_crit_avg_AIC_all(:,j) =  shuf_full_crit_avg.ModelCriterion.AIC;
% 
% %HEn1 ME
% shuf_spec_full_he_me = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  shuf_HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_he_me = fitlme(T_Both,shuf_spec_full_he_me);
% shuf_full_he_me_coef_all(:,j) = shuf_full_he_me.Coefficients.Estimate;
% shuf_full_he_me_AIC_all(:,j) =  shuf_full_he_me.ModelCriterion.AIC;
% %Hen1:n1
% shuf_spec_full_he_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:shuf_HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_he_n1 = fitlme(T_Both,shuf_spec_full_he_n1);
% shuf_full_he_n1_coef_all(:,j) = shuf_full_he_n1.Coefficients.Estimate;
% shuf_full_he_n1_AIC_all(:,j) =  shuf_full_he_n1.ModelCriterion.AIC;
% %hen1:MA
% shuf_spec_full_he_avg = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:shuf_HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% shuf_full_he_avg = fitlme(T_Both,shuf_spec_full_he_avg);
% shuf_full_he_avg_coef_all(:,j) = shuf_full_he_avg.Coefficients.Estimate;
% shuf_full_he_avg_AIC_all(:,j) =  shuf_full_he_avg.ModelCriterion.AIC;

end
%% shuffle v acutal permutation tests for simple lme
%asking, how many times is there a value in the shuffled daataset that is
%at least as extreme as the actually obtained value? Set alpha to .05, so
%if, for instance, we have 1000 shuffled datasets, we say there can be no
%more than 50 times that a shuffled variable is larger than the actual.
n_p_n1_simple = sum( abs(shuf_n1_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(2)))/size(shuf_n1_simple_coef_all,2);
n_p_n2_simple = sum( abs(shuf_n2_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(3)))/size(shuf_n2_simple_coef_all,2);
n_p_n3_simple = sum( abs(shuf_n3_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(4)))/size(shuf_n3_simple_coef_all,2);
n_p_n4_simple = sum( abs(shuf_n4_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(5)))/size(shuf_n4_simple_coef_all,2);
n_p_n5_simple = sum( abs(shuf_n5_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(6)))/size(shuf_n5_simple_coef_all,2);
n_p_n6_simple = sum( abs(shuf_n6_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(7)))/size(shuf_n6_simple_coef_all,2);
n_p_n7_simple = sum( abs(shuf_n7_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(8)))/size(shuf_n7_simple_coef_all,2);
n_p_n8_simple = sum( abs(shuf_n8_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(9)))/size(shuf_n8_simple_coef_all,2);
n_p_n9_simple = sum( abs(shuf_n9_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(10)))/size(shuf_n9_simple_coef_all,2);
n_p_n10_simple = sum( abs(shuf_n10_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(11)))/size(shuf_n10_simple_coef_all,2);
n_p_crit_simple = sum( abs(shuf_crit_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(12)))/size(shuf_crit_simple_coef_all,2);
n_p_ts_simple = sum( abs(shuf_ts_simple_coef_all(2,:)) > abs(lme_10_ma_totrew_coef(13)))/size(shuf_ts_simple_coef_all,2);
%put all the permutation pvals together in order of the coefs
n_p_1_to_10_totrew_simple = [n_p_n1_simple; n_p_n2_simple; n_p_n3_simple; n_p_n4_simple; n_p_n5_simple; n_p_n6_simple; n_p_n7_simple; n_p_n8_simple; n_p_n9_simple; n_p_n10_simple; n_p_crit_simple; n_p_ts_simple;];

%get means and sem for each nback
n1_shuf_mean = mean(shuf_n1_simple_coef_all(2,:)) ;
n1_shuf_sem = std(shuf_n1_simple_coef_all(2,:))/(sqrt(size( shuf_n1_simple_coef_all(2,:),2)));
n2_shuf_mean = mean(shuf_n2_simple_coef_all(2,:)) ;
n2_shuf_sem = std(shuf_n2_simple_coef_all(2,:))/(sqrt(size( shuf_n2_simple_coef_all(2,:),2)));
n3_shuf_mean = mean(shuf_n3_simple_coef_all(2,:)) ;
n3_shuf_sem = std(shuf_n3_simple_coef_all(2,:))/(sqrt(size( shuf_n3_simple_coef_all(2,:),2)));
n4_shuf_mean = mean(shuf_n4_simple_coef_all(2,:)) ;
n4_shuf_sem = std(shuf_n4_simple_coef_all(2,:))/(sqrt(size( shuf_n4_simple_coef_all(2,:),2)));
n5_shuf_mean = mean(shuf_n5_simple_coef_all(2,:)) ;
n5_shuf_sem = std(shuf_n5_simple_coef_all(2,:))/(sqrt(size( shuf_n5_simple_coef_all(2,:),2)));
n6_shuf_mean = mean(shuf_n6_simple_coef_all(2,:)) ;
n6_shuf_sem = std(shuf_n6_simple_coef_all(2,:))/(sqrt(size( shuf_n6_simple_coef_all(2,:),2)));
n7_shuf_mean = mean(shuf_n7_simple_coef_all(2,:)) ;
n7_shuf_sem = std(shuf_n7_simple_coef_all(2,:))/(sqrt(size( shuf_n7_simple_coef_all(2,:),2)));
n8_shuf_mean = mean(shuf_n8_simple_coef_all(2,:)) ;
n8_shuf_sem = std(shuf_n8_simple_coef_all(2,:))/(sqrt(size( shuf_n8_simple_coef_all(2,:),2)));
n9_shuf_mean = mean(shuf_n9_simple_coef_all(2,:)) ;
n9_shuf_sem = std(shuf_n9_simple_coef_all(2,:))/(sqrt(size( shuf_n9_simple_coef_all(2,:),2)));
n10_shuf_mean = mean(shuf_n10_simple_coef_all(2,:)) ;
n10_shuf_sem = std(shuf_n10_simple_coef_all(2,:))/(sqrt(size( shuf_n10_simple_coef_all(2,:),2)));
%put the mean and sem together in coef order for graping
n1_to_n10_shuf_mean = [n1_shuf_mean; n2_shuf_mean; n3_shuf_mean; n4_shuf_mean; n5_shuf_mean; n6_shuf_mean; n7_shuf_mean; n8_shuf_mean; n9_shuf_mean; n10_shuf_mean];
n1_to_n10_shuf_sem = [n1_shuf_sem; n2_shuf_sem; n3_shuf_sem; n4_shuf_sem; n5_shuf_sem; n6_shuf_sem; n7_shuf_sem; n8_shuf_sem; n9_shuf_sem; n10_shuf_sem; ];

%rather than the SEM of the shuffle coefs, lets calculaute the average SEM of the shuffled models 
n1_to_n10_shuf_avg_model_se = [mean(shuf_n1_simple_SE_all(2,:)); mean(shuf_n2_simple_SE_all(2,:)); mean(shuf_n3_simple_SE_all(2,:)); mean(shuf_n4_simple_SE_all(2,:)); mean(shuf_n5_simple_SE_all(2,:)); mean(shuf_n6_simple_SE_all(2,:)); mean(shuf_n7_simple_SE_all(2,:)); mean(shuf_n8_simple_SE_all(2,:)); mean(shuf_n9_simple_SE_all(2,:)); mean(shuf_n10_simple_SE_all(2,:)); ];
    
%% permutation tests of canonical model
%permutation tests from the AIC selected model
% n_p_full_n_perm_all= NaN(length(lme_reduced_coef),1);
% 
% %duration me shuffles
% n_p_full_n1_dur = sum( abs(shuf_full_n1dur_coef_all(2,:)) > abs(lme_reduced_coef(3)))/size(shuf_full_n1dur_coef_all,2);
% n_p_full_n_perm_all(3) = n_p_full_n1_dur;
% 
% n_p_full_n2_dur = sum( abs(shuf_full_n2dur_coef_all(2,:)) > abs(lme_reduced_coef(4)))/size(shuf_full_n2dur_coef_all,2);
% n_p_full_n_perm_all(4) = n_p_full_n2_dur;
% 
% n_p_full_n3_dur = sum( abs(shuf_full_n3dur_coef_all(2,:)) > abs(lme_reduced_coef(5)))/size(shuf_full_n3dur_coef_all,2);
% n_p_full_n_perm_all(5) = n_p_full_n3_dur;
% 
% n_p_full_n4_dur = sum( abs(shuf_full_n4dur_coef_all(2,:)) > abs(lme_reduced_coef(6)))/size(shuf_full_n4dur_coef_all,2);
% n_p_full_n_perm_all(6) = n_p_full_n4_dur;
% 
% n_p_full_n5_dur = sum( abs(shuf_full_n5dur_coef_all(2,:)) > abs(lme_reduced_coef(7)))/size(shuf_full_n5dur_coef_all,2);
% n_p_full_n_perm_all(7) = n_p_full_n5_dur;
% 
% n_p_full_n6_dur = sum( abs(shuf_full_n6dur_coef_all(2,:)) > abs(lme_reduced_coef(8)))/size(shuf_full_n6dur_coef_all,2);
% n_p_full_n_perm_all(8) = n_p_full_n6_dur;
% 
% n_p_full_ma_me = sum( abs(shuf_full_ma_coef_all(2,:)) > abs(lme_reduced_coef(14)))/size(shuf_full_ma_coef_all,2);
% n_p_full_n_perm_all(14) = n_p_full_ma_me;
% %also get mean and sem for the shuffled moving average for graphical purposes
% full_ma_me_mean = mean(shuf_full_ma_coef_all(2,:));
% full_ma_me_sem = std(shuf_full_ma_coef_all(2,:))/(sqrt(size( shuf_full_ma_coef_all(2,:),2)));
% 
% %now compare this to when we shuffle main effects and interactions
% %individually
% n_p_full_n_perm_all_indivshuffles =n_p_full_n_perm_all;
% 
% n_p_full_ipi_me = sum( abs(shuf_full_ipi_me_coef_all(2,:)) > abs(lme_reduced_coef(12)))/size(shuf_full_ipi_me_coef_all,2);
% n_p_full_ipi_int = sum( abs(shuf_full_ipi_n1_coef_all(15,:)) > abs(lme_reduced_coef(19)))/size(shuf_full_ipi_n1_coef_all,2);
% n_p_full_ipi_MA = sum( abs(shuf_full_ipi_avg_coef_all(21,:)) > abs(lme_reduced_coef(25)))/size(shuf_full_ipi_avg_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(12) = n_p_full_ipi_me;
% n_p_full_n_perm_all_indivshuffles(19) = n_p_full_ipi_int;
% n_p_full_n_perm_all_indivshuffles(25) = n_p_full_ipi_MA;
% 
% n_p_full_ipi2_me = sum( abs(shuf_full_ipi2_me_coef_all(2,:)) > abs(lme_reduced_coef(13)))/size(shuf_full_ipi2_me_coef_all,2);
% n_p_full_ipi2_int = sum( abs(shuf_full_ipi2_n2_coef_all(16,:)) > abs(lme_reduced_coef(20)))/size(shuf_full_ipi2_n2_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(13) = n_p_full_ipi2_me;
% n_p_full_n_perm_all_indivshuffles(20) = n_p_full_ipi2_int;
% 
% n_p_full_he_me = sum( abs(shuf_full_he_me_coef_all(2,:)) > abs(lme_reduced_coef(10)))/size(shuf_full_he_me_coef_all,2);
% n_p_full_he_int = sum( abs(shuf_full_he_n1_coef_all(15,:)) > abs(lme_reduced_coef(17)))/size(shuf_full_he_n1_coef_all,2);
% n_p_full_he_MA = sum( abs(shuf_full_he_avg_coef_all(21,:)) > abs(lme_reduced_coef(23)))/size(shuf_full_he_avg_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(10) = n_p_full_he_me;
% n_p_full_n_perm_all_indivshuffles(17) = n_p_full_he_int;
% n_p_full_n_perm_all_indivshuffles(23) = n_p_full_he_MA;
% 
% n_p_full_rew_me = sum( abs(shuf_full_rew_me_coef_all(2,:)) > abs(lme_reduced_coef(2)))/size(shuf_full_rew_me_coef_all,2);
% n_p_full_rew_int = sum( abs(shuf_full_rew_n1_coef_all(15,:)) > abs(lme_reduced_coef(15)))/size(shuf_full_rew_n1_coef_all,2);
% n_p_full_rew_MA = sum( abs(shuf_full_rew_avg_coef_all(21,:)) > abs(lme_reduced_coef(21)))/size(shuf_full_rew_avg_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(2) = n_p_full_rew_me;
% n_p_full_n_perm_all_indivshuffles(15) = n_p_full_rew_int;
% n_p_full_n_perm_all_indivshuffles(21) = n_p_full_rew_MA;
% 
% n_p_full_ts_me = sum( abs(shuf_full_ts_me_coef_all(2,:)) > abs(lme_reduced_coef(11)))/size(shuf_full_ts_me_coef_all,2);
% n_p_full_ts_int = sum( abs(shuf_full_ts_n1_coef_all(15,:)) > abs(lme_reduced_coef(18)))/size(shuf_full_ts_n1_coef_all,2);
% n_p_full_ts_MA = sum( abs(shuf_full_ts_avg_coef_all(21,:)) > abs(lme_reduced_coef(24)))/size(shuf_full_ts_avg_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(11) = n_p_full_ts_me;
% n_p_full_n_perm_all_indivshuffles(18) = n_p_full_ts_int;
% n_p_full_n_perm_all_indivshuffles(24) = n_p_full_ts_MA;
% 
% n_p_full_crit_me = sum( abs(shuf_full_crit_me_coef_all(2,:)) > abs(lme_reduced_coef(9)))/size(shuf_full_crit_me_coef_all,2);
% n_p_full_crit_int = sum( abs(shuf_full_crit_n1_coef_all(15,:)) > abs(lme_reduced_coef(16)))/size(shuf_full_crit_n1_coef_all,2);
% n_p_full_crit_MA = sum( abs(shuf_full_crit_avg_coef_all(21,:)) > abs(lme_reduced_coef(22)))/size(shuf_full_crit_avg_coef_all,2);
% n_p_full_n_perm_all_indivshuffles(9) = n_p_full_crit_me;
% n_p_full_n_perm_all_indivshuffles(16) = n_p_full_crit_int;
% n_p_full_n_perm_all_indivshuffles(22) = n_p_full_crit_MA;
% 
% %do it based on aic instead of coef value
% n_p_full_n_perm_all_AIC= NaN(length(lme_reduced_coef),1);
% lme_reduced_AIC = lme_reduced.ModelCriterion.AIC;
% 
% %with AIC, if the non-shuffled version is smaller, it is better. Thus we
% %want to see how many times the shuffled version of AIC is larger than the
% %acutal
% %if AIC in the shuffled is smaller - that means it's "better" 
% %that happens 12/100 times with n1. Thus, 12 of the shuffled ones we get a
% %better fit to teh data, and should be rejecteed
% n_p_full_n1_dur = sum( shuf_full_n1dur_AIC_all < lme_reduced_AIC)/size(shuf_full_n1dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(3) = n_p_full_n1_dur;
% 
% n_p_full_n2_dur = sum( abs(shuf_full_n2dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n2dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(4) = n_p_full_n2_dur;
% 
% n_p_full_n3_dur = sum( abs(shuf_full_n3dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n3dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(5) = n_p_full_n3_dur;
% 
% n_p_full_n4_dur = sum( abs(shuf_full_n4dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n4dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(6) = n_p_full_n4_dur;
% 
% n_p_full_n5_dur = sum( abs(shuf_full_n5dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n5dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(7) = n_p_full_n5_dur;
% 
% n_p_full_n6_dur = sum( abs(shuf_full_n6dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n6dur_AIC_all,2);
% n_p_full_n_perm_all_AIC(8) = n_p_full_n6_dur;
% 
% n_p_full_ma_me = sum( abs(shuf_full_ma_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ma_AIC_all,2);
% n_p_full_n_perm_all_AIC(14) = n_p_full_ma_me;
% 
% n_p_full_ipi_me = sum( abs(shuf_full_ipi_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_me_AIC_all,2);
% n_p_full_ipi_int = sum( abs(shuf_full_ipi_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_n1_AIC_all,2);
% n_p_full_ipi_MA = sum( abs(shuf_full_ipi_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_avg_AIC_all,2);
% n_p_full_n_perm_all_AIC(12) = n_p_full_ipi_me;
% n_p_full_n_perm_all_AIC(19) = n_p_full_ipi_int;
% n_p_full_n_perm_all_AIC(25) = n_p_full_ipi_MA;
% 
% n_p_full_ipi2_me = sum( abs(shuf_full_ipi2_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi2_me_AIC_all,2);
% n_p_full_ipi2_int = sum( abs(shuf_full_ipi2_n2_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi2_n2_AIC_all,2);
% n_p_full_n_perm_all_AIC(13) = n_p_full_ipi2_me;
% n_p_full_n_perm_all_AIC(20) = n_p_full_ipi2_int;
% 
% n_p_full_he_me = sum( abs(shuf_full_he_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_me_AIC_all,2);
% n_p_full_he_int = sum( abs(shuf_full_he_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_n1_AIC_all,2);
% n_p_full_he_MA = sum( abs(shuf_full_he_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_avg_AIC_all,2);
% n_p_full_n_perm_all_AIC(10) = n_p_full_he_me;
% n_p_full_n_perm_all_AIC(17) = n_p_full_he_int;
% n_p_full_n_perm_all_AIC(23) = n_p_full_he_MA;
% 
% n_p_full_rew_me = sum( abs(shuf_full_rew_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_me_AIC_all,2);
% n_p_full_rew_int = sum( abs(shuf_full_rew_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_n1_AIC_all,2);
% n_p_full_rew_MA = sum( abs(shuf_full_rew_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_avg_AIC_all,2);
% n_p_full_n_perm_all_AIC(2) = n_p_full_rew_me;
% n_p_full_n_perm_all_AIC(15) = n_p_full_rew_int;
% n_p_full_n_perm_all_AIC(21) = n_p_full_rew_MA;
% 
% n_p_full_ts_me = sum( abs(shuf_full_ts_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_me_AIC_all,2);
% n_p_full_ts_int = sum( abs(shuf_full_ts_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_n1_AIC_all,2);
% n_p_full_ts_MA = sum( abs(shuf_full_ts_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_avg_AIC_all,2);
% n_p_full_n_perm_all_AIC(11) = n_p_full_ts_me;
% n_p_full_n_perm_all_AIC(18) = n_p_full_ts_int;
% n_p_full_n_perm_all_AIC(24) = n_p_full_ts_MA;
% 
% n_p_full_crit_me = sum( abs(shuf_full_crit_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_me_AIC_all,2);
% n_p_full_crit_int = sum( abs(shuf_full_crit_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_n1_AIC_all,2);
% n_p_full_crit_MA = sum( abs(shuf_full_crit_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_avg_AIC_all,2);
% n_p_full_n_perm_all_AIC(9) = n_p_full_crit_me;
% n_p_full_n_perm_all_AIC(16) = n_p_full_crit_int;
% n_p_full_n_perm_all_AIC(22) = n_p_full_crit_MA;
%      

%% Predictions 

%predict LP_Durations from the model (predicted response using model data by default)
 [yhat yhatCI yhatDF]= predict(lme_reduced,'Simultaneous',true);
%   [yhat yhatCI yhatDF]= predict(lme_reduced,'Simultaneous',false);

%Model that is LM instead of LME
%  [yhat yhatCI]= predict(lm_reduced,T10_Logical_and_Continuous,'Simultaneous',true);

%calculate how often the actual value lies between the CI for the
 %predicted
 correctish_pred =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI(:,1);
 correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.LP_Durations_All);
%24.1% in 33-3 beh
%LM predicts 11.8%
%LM trained on 33-3 predicts 11.6% of 32-2 

%predict a smoothed dataset
%smooth the 10 presses to the left
 smootheddata = smoothdata(T10_Logical_and_Continuous.LP_Durations_All,'gaussian',[0 10]);
 correctish_pred_smooth =smootheddata<=yhatCI(:,2) & smootheddata >=yhatCI(:,1);
 correctCI_prop_Smooth = sum(correctish_pred_smooth)/length(smootheddata) ;
 %41.2% in 33-3 beh
 
%What about predicting order shuffled Durations?
 correctish_shuf_pred =Shuffled_Durations_All <=yhatCI(:,2) & Shuffled_Durations_All >=yhatCI(:,1);
 correctCI_shuf_prop = sum(correctish_shuf_pred)/length(Shuffled_Durations_All);
 %22.8% in 33-3 beh
 
%so it turns out it predicts about the same number of presses (23% in shuf
%vs 25% in actual. This likely reflects the fact that the prediction conf
%bounds dont change drastically across the session. Now check to see how
%often THE SAME press is identified in the shuffle and actual
actual_is_shuf = correctish_shuf_pred == 1 & correctish_pred ==1;
sum(actual_is_shuf)/length(actual_is_shuf);
 %Much Smaller: 6.7% are still predicted. Thus, we see a difference of
 %22.55 - 5.9 = 16.1% more presses predicted
 
 %How well does the model simply predict an increase a decrease relative to
%the preceding press?
diff_yhat = diff(yhat);
diff_x = diff(T10_Logical_and_Continuous.LP_Durations_All);
dixx_x_yhat = [diff_x diff_yhat];
%convert the numbers to either -1 or +1 to see how often it predicts a
%increase or decrease accurately
 dixx_x_yhat(dixx_x_yhat > 0) =1;
 dixx_x_yhat(dixx_x_yhat < 0) =-1;

inc_pred =0;
dec_pred = 0;
for i = 1:length(dixx_x_yhat)
    if dixx_x_yhat(i,1) == 1 && dixx_x_yhat(i,2) == 1
        inc_pred = inc_pred+1;
    elseif  dixx_x_yhat(i,1) == -1 && dixx_x_yhat(i,2) == -1
        dec_pred = dec_pred +1;
    end
    
end

actual_inc = sum(dixx_x_yhat(:,1)==1);
correct_inc_pred = inc_pred/actual_inc;
 actual_dec=sum(dixx_x_yhat(:,1)==-1);
correct_dec_pred =dec_pred/actual_dec;

sum(dixx_x_yhat(:,1) == dixx_x_yhat(:,2))/length(dixx_x_yhat);
%Predicts 43.1% of the time when presses increase/decrease overall -
%roughly the same for both increasing and decreasing individually

%How well does the model accurately predict rewarded/unrewarded presses?
%use the upper CI to see if the prediciton would be rewarded or not
%need numerical crteria% indicator for comparisons of duration to criterion
T10_Logical_and_Continuous.criteria_indicator_All_numerical = criteria_indicator_All;
pred_rew = yhatCI(:,2) >= T10_Logical_and_Continuous.criteria_indicator_All_numerical; 
pred_rew = double(pred_rew);
% ll here is the n-0 reward indicator var
rew_red_actual = [logical_LP_n_back_DataRaster(:,11) pred_rew];

success_pred =0;
fail_pred = 0;
for i = 1:length(rew_red_actual)
%     for i = (41267:414460 )
%if both the actual and predicted reward are 1
    if rew_red_actual(i,1) == 1 && rew_red_actual(i,2) == 1
        %then iterate the accurate success_pred variable
        success_pred = success_pred+1;
        %ditto for if it agrees on a failed press
    elseif  rew_red_actual(i,1) == 0 && rew_red_actual(i,2) == 0
        fail_pred = fail_pred +1;
    end
    
end
actual_rewards = sum(logical_LP_n_back_DataRaster(:,11));
correct_reward_pred = success_pred/actual_rewards;
%47.3% of rewards predicted accurately

% correct_reward_pred = success_pred/187;

 fail_idxs=sum(logical_LP_n_back_DataRaster(:,11) == 0);
correct_fail_pred =fail_pred/fail_idxs;
%81.3% of failed presses predicted accurately
% correct_fail_pred =fail_pred/187

sum(pred_rew == logical_LP_n_back_DataRaster(:,11))/length(pred_rew);
%73.8% of presses predicted accurately in terms of success/failure
%lM predicts 75.7% in 33-3
%simple lm trained on 33-3 predicts 78.1% of 32-2
%% predict with the simple n10 and control variable model
 [yhat_simple yhatCI_simple yhatDF_simple] = predict(lme_10_ma_totrew,'Simultaneous',true);
%        
%calculate how often the actual value lies between the CI for the
 %predicted
 correctish_pred_simple =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI_simple(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI_simple(:,1);
 correctCI_prop_simple = sum(correctish_pred_simple)/length(T10_Logical_and_Continuous.LP_Durations_All);
 %17.7% of presses predicted accurately by the simple model, which is worse
 %than the more complex model performacne of 22.5%. Good that the more
 %complex model performs better
 

%% Stimulation and probabilistic LMEs
%comment out unless looking at opto data or probabalistic reward

%% Stimulation during the lever press
lme3_time_stim_n_int_ctl_ma_int_full = fitlme(T10_Logical_and_Continuous,'LP_Durations_All ~stim_indicator_All*n_minus_one_Durations_All + stim_indicator_n1_All*n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All + criteria_percent_indicator+ipi2*n_minus_two_Durations_All+n_minus_one_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n_1_Indicator_All+LP_Timestamps_All*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All + ipi1:moving_average_lp_length_n7andback_All +  moving_average_lp_length_n7andback_All:n_minus_one_All+moving_average_lp_length_n7andback_All:criteria_percent_indicator  + LP_Timestamps_All:moving_average_lp_length_n7andback_All+stim_indicator_All*moving_average_lp_length_n7andback_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All  + (1|mouse_indicator)+(1|day_indicator)');
lme3_time_stim_n_int_ctl_ma_int_full_coef = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Estimate;
lme3_time_stim_n_int_ctl_ma_int_full_se = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.SE;
lme3_time_stim_n_int_ctl_ma_int_full_pval = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.pValue;
lme3_time_stim_n_int_ctl_ma_int_full_names = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Name;
lme3_time_stim_n_int_ctl_ma_int_full_df = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.DF;
lme3_time_stim_n_int_ctl_ma_int_full_tstat = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.tStat;
lme3_time_stim_n_int_ctl_ma_int_full_upper = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Upper;
lme3_time_stim_n_int_ctl_ma_int_full_lower = lme3_time_stim_n_int_ctl_ma_int_full.Coefficients.Lower;
lme3_time_stim_n_int_ctl_ma_int_full_anova = anova(lme3_time_stim_n_int_ctl_ma_int_full);
lme3_time_stim_n_int_ctl_ma_int_full_fstat = lme3_time_stim_n_int_ctl_ma_int_full_anova.FStat;
lme3_time_stim_n_int_ctl_ma_int_full_fpval = lme3_time_stim_n_int_ctl_ma_int_full_anova.pValue;

%% stimulation after lmes
model_spec_reduced_stim = 'LP_Durations_All ~ stim_indicator_n1_All*n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All + criteria_percent_indicator+ipi2*n_minus_two_Durations_All+n_minus_one_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n_1_Indicator_All+LP_Timestamps_All*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All + ipi1:moving_average_lp_length_n7andback_All +  moving_average_lp_length_n7andback_All:n_minus_one_All+moving_average_lp_length_n7andback_All:criteria_percent_indicator  + LP_Timestamps_All:moving_average_lp_length_n7andback_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All  + (1|mouse_indicator)+(1|day_indicator)';
lme_reduced_stim = fitlme(T10_Logical_and_Continuous,model_spec_reduced_stim)
lme_reduced_stim_name = lme_reduced_stim.Coefficients.Name;
lme_reduced_stim_coef = lme_reduced_stim.Coefficients.Estimate;
lme_reduced_stim_SE = lme_reduced_stim.Coefficients.SE;
lme_reduced_stim_pval = lme_reduced_stim.Coefficients.pValue;
lme_reduced_stim_anova = anova(lme_reduced_stim);
lme_reduced_stim_fstat = lme_reduced_stim_anova.FStat;
lme_reduced_stim_fpval = lme_reduced_stim_anova.pValue;

 %% for probability add in main effect and interaction of success + rew vs success no reward
% remove n_minus_one_all effects, since it only cares whether duration was long enough, not whether there was reward
% (this only matters if reward is probbalistic)
% instead use the reward_indicator_n1_All variable
% where 0 = fail, no reward, 1 = success, no reward, and 2 = success, yes
% reward
T10_Logical_and_Continuous.reward_indicator_n1_All = categorical(T10_Logical_and_Continuous.reward_indicator_n1_All);
T10_Logical_and_Continuous.reward_indicator_All = categorical(T10_Logical_and_Continuous.reward_indicator_All);

%BIC selected model
model_spec_moderate_remove_proability = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

lme__moderate_remove_probability = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability)
lme__moderate_remove_probability_coef = lme__moderate_remove_probability.Coefficients.Estimate;
lme__moderate_remove_probability_pvals = lme__moderate_remove_probability.Coefficients.pValue;	
lme__moderate_remove_probability_se = lme__moderate_remove_probability.Coefficients.SE;
lme__moderate_remove_probability_names = lme__moderate_remove_probability.Coefficients.Name;
prob_anova = anova(lme__moderate_remove_probability);
prob_anova_F = prob_anova.FStat;
prob_anova_pval = prob_anova.pValue;
[b bnames bstats] = randomEffects(lme__moderate_remove_probability);


% For running all the prob groups together, we can tag the animals that
% were exposed to 25, 50, and 75% reward respectively
%The mouse_indicator variable just runs through all files in numerical
%order, so need to convert experimental mouse numbers into matlab mouse
%indicator. Basically, no mouse 6, so need to remove 1 from all after 6
%matlab number : mouse number
% 1 = 1, 2 = 2, 3 = 3, 4 = 4, 5 = 5, 6 = 7 ,7 = 8 , 8 = 9, 9 = 10, 10 = 11,
% 11 = 12, 12 = 13 13 = 14, 14 = 15, 15 = 16
%25% = 1,3,4,8,10
%50% = 2,7,9,14,15
%75% = 5,6,11,12,13
%25%
prob_indicator = nan(length(T10_Logical_and_Continuous.mouse_indicator ),1);
indicator_25 = T10_Logical_and_Continuous.mouse_indicator == '1'|T10_Logical_and_Continuous.mouse_indicator == '3'| T10_Logical_and_Continuous.mouse_indicator == '4' |T10_Logical_and_Continuous.mouse_indicator == '8'| T10_Logical_and_Continuous.mouse_indicator == '10';
indicator_50 = T10_Logical_and_Continuous.mouse_indicator == '2'|T10_Logical_and_Continuous.mouse_indicator == '7'| T10_Logical_and_Continuous.mouse_indicator == '9' |T10_Logical_and_Continuous.mouse_indicator == '14'| T10_Logical_and_Continuous.mouse_indicator == '15';
indicator_75 = T10_Logical_and_Continuous.mouse_indicator == '5'|T10_Logical_and_Continuous.mouse_indicator == '6'| T10_Logical_and_Continuous.mouse_indicator == '11' |T10_Logical_and_Continuous.mouse_indicator == '12'| T10_Logical_and_Continuous.mouse_indicator == '13';
prob_indicator(indicator_25) = 25;
prob_indicator(indicator_50) = 50;
prob_indicator(indicator_75) = 75;
T10_Logical_and_Continuous.prob_indicator = categorical(prob_indicator);

sum(T10_Logical_and_Continuous.prob_indicator =='75');

%BIC selected model
%include interactions with prob group and reward indicator, and 3-way with
%prob_group:reward:n1 dur
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

lme__moderate_remove_probability_probgroup = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability_probgroup,'DummyVarCoding','reference')
lme__moderate_remove_probability_probgroup_coef = lme__moderate_remove_probability_probgroup.Coefficients.Estimate;
lme__moderate_remove_probability_probgroup_pvals = lme__moderate_remove_probability_probgroup.Coefficients.pValue;	
lme__moderate_remove_probability_probgroup_se = lme__moderate_remove_probability_probgroup.Coefficients.SE;
lme__moderate_remove_probability_probgroup_names = lme__moderate_remove_probability_probgroup.Coefficients.Name;
[b bnames bstats] = randomEffects(lme__moderate_remove_probability_probgroup);

probgroups_anova = anova(lme__moderate_remove_probability_probgroup);
probgroups_anova_F = probgroups_anova.FStat;
probgroups_anova_pval = probgroups_anova.pValue;

%% Save all the data for later use if wanted. Need to save it in old matlab format
%file format (v7.3) due to size.
% save('1000shufdata-33-3' , '-v7.3')

toc   