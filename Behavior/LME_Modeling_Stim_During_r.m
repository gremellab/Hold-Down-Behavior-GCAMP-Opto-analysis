%% Purpose
%The main purpose of this script is to create Linear Mixed Effect Models (LMEs)
%to predict the duration of lever presses (n) given subjective experiential
%information. This includes things like n-back press durations, n-back
%headentries (HE), n-back reward, interpressinterval (IPI), time in
%session, % of presses over criteria, cumulative rewards across a session.
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
%utilize the Data structure obtained from the MEDPC_Behavior_Extract_For_Regressions
%This includes things like generating within session performance graphs,
%histograms of press durations, press duration entropy, and raw duration
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

load('67-1-67-2-opto-800ms1-6') %load the data structure
% total_mice =16; %hard code number of mice if needed, otherwise get it via
% the size of the data structure

%For the behavioral dataset, we had two mice that made either 0 or 1 LP on
%day 1. We need to figure out how to handle uneven mice across days

day1_total_mice = size(Data.Day(1).Mouse,2);
total_mice = size(Data.Day(2).Mouse,2);
%get number of days
day_number = size(Data.Day,2);
 
%Initialize Variables
LP_Durations_Cell = {};
LP_Timestamps_Cell ={};
Rewarded_Durations_Cell = {};
Unrewarded_Durations_Cell ={};
IPI_Cell ={};
Total_HE_Cell =[];
Shuffled_LP_Duration_Cell ={};
Total_HE_Prop_Cell =[];
Rewarded_HE_Prop_Cell = [];
Unrewarded_HE_Prop_Cell =[];
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
total_reward_indicator_Cell ={};
total_reward_indicator_prob_Cell ={};
reward_indicator_Cell ={};
reward_indicator_n1_Cell ={};
reward_indicator_n2_Cell ={};
reward_indicator_n3_Cell ={};
reward_indicator_n4_Cell ={};
moving_average_lp_length_n7andback_Cell ={};
lp_ent_Cell =[];
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
 
 CPD_n1_criteria_simple_Cell ={};
 CPD_n1_criteria_complex_Cell ={};
 CPD_MA_criteria_complex_Cell ={};
 CPD_n1rew_criteria_complex_Cell ={};
 total_time_and_rewards_Cell ={};
 
 mean_Shuffled_Length_After_Success_All_Shuffles =[];
 mean_Shuffled_Length_After_Failure_All_Shuffles =[];
 std_Shuffled_Length_After_Success_All_Shuffles =[];
 std_Shuffled_Length_After_Failure_All_Shuffles =[];
 mean_Shuffled_Change_After_Success_All_Shuffles =[];
 mean_Shuffled_Change_After_Failure_All_Shuffles =[];
%  entropy_Shuffled_Change_After_Success_All_Shuffles =[];
%  entropy_Shuffled_Change_After_Failure_All_Shuffles =[];

%% Loop through data structure to get data input into Cells
for i = 1:day_number %outer loop goes across days
%     if i == 1 %on day 1, we had 2 mice with only 1 lp, exclude these mice from day 1
%         total_mice = day1_total_mice;
%     else
%         total_mice = size(Data.Day(2).Mouse,2);
%     end

    for j = 1:total_mice %inner loop goes across mice within a day
         mean_Shuffled_Length_After_Success_All_Shuffles = [mean_Shuffled_Length_After_Success_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Length_After_Success_Dist];
         mean_Shuffled_Length_After_Failure_All_Shuffles = [mean_Shuffled_Length_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Length_After_Failure_Dist];
         std_Shuffled_Length_After_Success_All_Shuffles = [std_Shuffled_Length_After_Success_All_Shuffles; Data.Day(i).Mouse(j).std_Shuffled_Lever_Press_Length_After_Success_Dist];
         std_Shuffled_Length_After_Failure_All_Shuffles = [std_Shuffled_Length_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).std_Shuffled_Lever_Press_Length_After_Failure_Dist];
   
         mean_Shuffled_Change_After_Success_All_Shuffles = [mean_Shuffled_Change_After_Success_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Change_After_Success_Dist];
         mean_Shuffled_Change_After_Failure_All_Shuffles = [mean_Shuffled_Change_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).mean_Shuffled_Lever_Press_Change_After_Failure_Dist]; 
%          entropy_Shuffled_Change_After_Success_All_Shuffles =[entropy_Shuffled_Change_After_Success_All_Shuffles; Data.Day(i).Mouse(j).entropy_Shuffled_Lever_Press_Length_After_Success_Dist];
%          entropy_Shuffled_Change_After_Failure_All_Shuffles =[entropy_Shuffled_Change_After_Failure_All_Shuffles; Data.Day(i).Mouse(j).entropy_Shuffled_Lever_Press_Length_After_Failure_Dist];
         total_time_and_rewards_Cell =[total_time_and_rewards_Cell Data.Day(i).Mouse(j).total_time_and_rewards];
        
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
     CPD_n1_criteria_simple_Cell =[CPD_n1_criteria_simple_Cell Data.Day(i).Mouse(j).indiv_CPD_n1_criteria_simple];
     CPD_n1_criteria_complex_Cell =[CPD_n1_criteria_complex_Cell Data.Day(i).Mouse(j).indiv_CPD_n1_criteria_complex ];
     CPD_MA_criteria_complex_Cell =[CPD_MA_criteria_complex_Cell Data.Day(i).Mouse(j).indiv_CPD_MA_criteria_complex ];
     CPD_n1rew_criteria_complex_Cell = [CPD_n1rew_criteria_complex_Cell  Data.Day(i).Mouse(j).indiv_CPD_n1rew_criteria_complex];
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
     total_reward_indicator_Cell = [total_reward_indicator_Cell  Data.Day(i).Mouse(j).total_reward_indicator]; 
     total_reward_indicator_prob_Cell = [total_reward_indicator_prob_Cell  Data.Day(i).Mouse(j).total_reward_indicator_prob]; 
     lp_ent_Cell = [lp_ent_Cell Data.Day(i).Mouse(j).lpent];
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
     IPI_Cell = [IPI_Cell Data.Day(i).Mouse(j).Lever_Press_IPI];
     Total_HE_Cell = [Total_HE_Cell Data.Day(i).Mouse(j).Total_Head_Entries];
     Shuffled_LP_Duration_Cell =[Shuffled_LP_Duration_Cell Data.Day(i).Mouse(j).Shuffled_Lever_Press_lengths];
     Total_HE_Prop_Cell =[Total_HE_Prop_Cell Data.Day(i).Mouse(j).Total_HE_Prop];
     Rewarded_HE_Prop_Cell = [Rewarded_HE_Prop_Cell Data.Day(i).Mouse(j).Reward_HE_Prop];
     Unrewarded_HE_Prop_Cell =[Unrewarded_HE_Prop_Cell Data.Day(i).Mouse(j).Unreward_HE_Prop];
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
              
%reshape the cells to get them into cells that are num mice(row) by num day
%(col). This is used primarilyfor graphing purposes
Total_HE_Prop_By_Mouse = reshape(Total_HE_Prop_Cell,[total_mice,day_number]);
Reward_HE_Prop_By_Mouse = reshape(Rewarded_HE_Prop_Cell,[total_mice,day_number]);
Unrewarded_HE_Prop_By_Mouse = reshape(Unrewarded_HE_Prop_Cell,[total_mice,day_number]);
Total_HE_Prop_By_Mouse_mean = nanmean(Total_HE_Prop_By_Mouse,2);
Reward_HE_Prop_By_Mouse_mean = nanmean(Reward_HE_Prop_By_Mouse,2);
Unrewarded_HE_Prop_By_Mouse_mean = nanmean(Unrewarded_HE_Prop_By_Mouse,2);

IPI_by_mouse =  reshape(IPI_Cell,[total_mice,day_number]);
IPI_by_mouse = cellfun(@nanmean,IPI_by_mouse);
IPI_by_mouse_med =  reshape(IPI_Cell,[total_mice,day_number]);
IPI_by_mouse_med = cellfun(@nanmedian,IPI_by_mouse_med);
lp_ent_Cell_byday =reshape(lp_ent_Cell,[total_mice,day_number]);


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
% 
% after_failure_duration_entropy = cellfun(@quickEnt,after_failure_duration_By_Mouse);
% after_failure_duration_entropy_overallmouse = mean(after_failure_duration_entropy,2);
% after_success_duration_entropy =  cellfun(@quickEnt,after_success_duration_By_Mouse);
% after_success_duration_entropy_overallmouse = mean(after_success_duration_entropy,2);
% success_failure_entropy_diff = after_success_duration_entropy - after_failure_duration_entropy;
% success_failure_entropy_diff_overallmouse = nanmean(success_failure_entropy_diff,2);
% success_failure_entropy_diff_negative =success_failure_entropy_diff < 0;
% prop_neg_success_failure_entropy_diff = sum(success_failure_entropy_diff_negative,2)/length(success_failure_entropy_diff_negative);

%shuffled after success entropy
% average_shuffled_after_success_entropy = mean(entropy_Shuffled_Change_After_Success_All_Shuffles,2);
% average_shuffled_after_success_entropy_byday = reshape(average_shuffled_after_success_entropy,[total_mice,day_number]);
% average_shuffled_after_success_entropy_overallmouse = mean(average_shuffled_after_success_entropy_byday,2);
% 
%how often is there less entropy in the actual vs the shuffled
% success_entropy_v_shuffled_pval=[];
% after_success_duration_entropy_inline = reshape(after_success_duration_entropy,[size(after_success_duration_entropy,1)*size(after_success_duration_entropy,2),1]);
% for perm_it = 1:length(after_success_duration_entropy_inline)
%     success_entropy_v_shuffled_pval =[success_entropy_v_shuffled_pval; sum(average_shuffled_after_success_entropy(perm_it,:) <  after_success_duration_entropy_inline(perm_it))/size(average_shuffled_after_success_entropy,2)];
% end
% success_entropy_v_shuffled_pval_bymouseday = reshape(success_entropy_v_shuffled_pval,[total_mice,day_number]);
% sum(success_entropy_v_shuffled_pval_bymouseday <= .05,2)/day_number


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
sum(success_std_v_shuffled_pval_bymouseday <= .05,2)/day_number


 overall_shuffled_after_failure_mean = mean(mean_Shuffled_Length_After_Failure_All_Shuffles,2);
 overall_shuffled_after_success_mean_by_mouse =reshape(overall_shuffled_after_success_mean,[total_mice,day_number]);
 overall_shuffled_after_failure_mean_by_mouse =reshape(overall_shuffled_after_failure_mean,[total_mice,day_number]);
 %Now we want to see how many times the mean in actual is greater than that
%in the shuffled
success_mean_v_shuffled_pval=[];
after_success_duration_mean_inline = reshape(after_success_duration_mean,[size(after_success_duration_mean,1)*size(after_success_duration_mean,2),1]);
for perm_it = 1:length(after_success_duration_std_inline)
    success_mean_v_shuffled_pval =[success_mean_v_shuffled_pval; sum(overall_shuffled_after_success_mean(perm_it,:) <  after_success_duration_mean_inline(perm_it))/size(after_success_duration_mean_inline,2)];
end
success_mean_v_shuffled_pval_bymouseday = reshape(success_mean_v_shuffled_pval,[total_mice,day_number]);
sum(success_mean_v_shuffled_pval_bymouseday <= .05,2)/day_number

 
 
 
 average_shuffled_after_failure_std = mean(std_Shuffled_Length_After_Failure_All_Shuffles,2); 
 average_shuffled_after_failure_std_byday = reshape(average_shuffled_after_failure_std ,[total_mice,day_number]);
average_shuffled_after_failure_std_bymouse =mean(average_shuffled_after_failure_std_byday,2);
 
 
success_failure_shuffled_diff =  overall_shuffled_after_success_mean_by_mouse-overall_shuffled_after_failure_mean_by_mouse;
success_failure_shuffled_diff_overallmousemean = mean(success_failure_shuffled_diff,2); 
%one-sided permutation test to see how many times there is a value as
%extremely negative in actual vs shuffled success/failure differences. This
%is to look for mice that make short LPs after reward
%calculate the mean diff per animal/day
after_failure_success_diff_Cell = cellfun(@mean,after_success_duration_Cell)' - cellfun(@mean,after_failure_duration_Cell)';
success_fail_diff_pval=[];
for perm_it = 1:length(after_failure_success_diff_Cell)
    success_fail_diff_pval =[success_fail_diff_pval; sum(overall_shuffled_after_success_failure_diff(perm_it,:) <  after_failure_success_diff_Cell(perm_it))/size(overall_shuffled_after_success_failure_diff,2)];
end

success_fail_diff_pval_by_mouseday = reshape(success_fail_diff_pval,[total_mice,day_number]);

%Find the number of days out of the total where a mouse fails this check
sum(success_fail_diff_pval_by_mouseday <= .05,2)/day_number
%the 5th mouse in 33-3 (numbered as mouse 7 in the datafiles) exceeds this
%criterion on 50% of days

% %Now lets look at change after success/failure
%  mean_Shuffled_Change_After_Success_All_Shuffles =[];
%  mean_Shuffled_Change_After_Failure_All_Shuffles =[];
 after_failure_change_mean = cellfun(@nanmean,after_failure_change_Cell);
 after_failure_change_mean = after_failure_change_mean';
  after_failure_change_mean_bymouseday = reshape(after_failure_change_mean,[total_mice,day_number]);
  after_failure_change_mean_mouse_overall = mean(after_failure_change_mean_bymouseday,2);
  after_success_change_mean = cellfun(@nanmean,after_success_change_Cell);
 after_success_change_mean = after_success_change_mean';
  after_success_change_mean_bymouseday = reshape(after_success_change_mean,[total_mice,day_number]);
  after_success_change_mean_mouse_overall = mean(after_success_change_mean_bymouseday,2);
  
 overall_mean_Shuffled_Change_After_Success_All_Shuffles = mean(mean_Shuffled_Change_After_Success_All_Shuffles,2); 
 overall_mean_Shuffled_Change_After_Failure_All_Shuffles = mean(mean_Shuffled_Change_After_Failure_All_Shuffles,2); 
 mean_of_success_change_shuffles_byday = reshape(overall_mean_Shuffled_Change_After_Success_All_Shuffles,[total_mice,day_number]);
 mean_of_failure_change_shuffles_byday = reshape(overall_mean_Shuffled_Change_After_Failure_All_Shuffles,[total_mice,day_number]); 
 
change_after_success_shuffled_mouse_mean = mean(mean_of_success_change_shuffles_byday,2);
 change_after_failure_shuffled_mouse_mean = mean(mean_of_failure_change_shuffles_byday,2);
 
actual_shuffled_diff_success_change_byday = abs( after_success_change_mean_bymouseday) - abs(mean_of_success_change_shuffles_byday);
 
 actual_shuffled_diff_success_change = abs(after_success_change_mean_mouse_overall) - abs(change_after_success_shuffled_mouse_mean);
 actual_shuffled_diff_failure_change = after_failure_change_mean_mouse_overall - change_after_failure_shuffled_mouse_mean;
 
 %check how many days per mouse they had a larger actual success change
 %than as predicted from the shuffled distributions
 success_change_actual_v_shuffled = mean_of_success_change_shuffles_byday > after_success_change_mean_bymouseday;
prop_neg_success_actual_v_shuffled_diff = sum(success_change_actual_v_shuffled,2)/length(success_change_actual_v_shuffled);

 
 %how often is the change in duration after success greater in actual v
 %shuffled
success_change_pval=[];
for perm_it = 1:length(after_success_change_mean)
    success_change_pval =[success_change_pval; sum(mean_Shuffled_Change_After_Success_All_Shuffles(perm_it,:) <  after_success_change_mean(perm_it))/size(overall_shuffled_after_success_failure_diff,2)];
end
success_change_pval_by_mouseday = reshape(success_change_pval,[total_mice,day_number]);
sum(success_change_pval_by_mouseday <= .05,2)/day_number

%how often is change in duration after failure greater in actual v shuffled
fail_change_pval=[];
for perm_it = 1:length(after_failure_change_mean)
    fail_change_pval =[fail_change_pval; sum(mean_Shuffled_Change_After_Failure_All_Shuffles(perm_it,:) >  after_failure_change_mean(perm_it))/size(overall_shuffled_after_success_failure_diff,2)];
end
fail_change_pval_by_mouseday = reshape(fail_change_pval,[total_mice,day_number]);
sum(fail_change_pval_by_mouseday <= .05,2)/day_number



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
total_reward_indicator_All=[];
total_reward_indicator_prob_All=[];

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

 CPD_n1_criteria_simple_All =[];
 CPD_n1_criteria_complex_All = [];
 CPD_MA_criteria_complex_All = [];
 CPD_n1rew_criteria_complex_All = [];
 total_time_and_rewards_All = [];

for i = 1:length(LP_Durations_Cell)  
total_time_and_rewards_All =[total_time_and_rewards_All; total_time_and_rewards_Cell{i}];
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
 CPD_n1_criteria_simple_All =[CPD_n1_criteria_simple_All; CPD_n1_criteria_simple_Cell{i}];
 CPD_n1_criteria_complex_All =[CPD_n1_criteria_complex_All; CPD_n1_criteria_complex_Cell{i} ];
 CPD_MA_criteria_complex_All =[CPD_MA_criteria_complex_All; CPD_MA_criteria_complex_Cell{i}];
 CPD_n1rew_criteria_complex_All = [CPD_n1rew_criteria_complex_All; CPD_n1rew_criteria_complex_Cell{i}];
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
     total_reward_indicator_All = [total_reward_indicator_All; total_reward_indicator_Cell{i}];
          total_reward_indicator_prob_All = [total_reward_indicator_prob_All; total_reward_indicator_prob_Cell{i}];

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
T10_Logical_and_Continuous.ipi7 = ipi_n_back_matrix(:,7);

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
T10_Logical_and_Continuous.moving_mean=moving_mean_All;
T10_Logical_and_Continuous.moving_average_lp_length_n7andback_All=moving_average_lp_length_n7andback_All;
%LP_All is logical for reward, and then n-back logicals
T10_Logical_and_Continuous.LP_All =categorical(T10_Logical_and_Continuous.LP_All);
T10_Logical_and_Continuous.n_minus_one_All =categorical(T10_Logical_and_Continuous.n_minus_one_All);
T10_Logical_and_Continuous.n_minus_two_All =categorical(T10_Logical_and_Continuous.n_minus_two_All);
T10_Logical_and_Continuous.n_minus_three_All =categorical(T10_Logical_and_Continuous.n_minus_three_All);
T10_Logical_and_Continuous.n_minus_four_All =categorical(T10_Logical_and_Continuous.n_minus_four_All);
T10_Logical_and_Continuous.n_minus_five_All =categorical(T10_Logical_and_Continuous.n_minus_five_All);
T10_Logical_and_Continuous.n_minus_six_All =categorical(T10_Logical_and_Continuous.n_minus_six_All);

%create an indicator for the last 3 days of training as a "late idx" if
%interested in seeing how late in training affects the model
T10_Logical_and_Continuous.day_indicator = double(T10_Logical_and_Continuous.day_indicator);
lateidx=T10_Logical_and_Continuous.day_indicator == 12:14;
lateidx=any(lateidx,2);
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.lateidx = lateidx;
%indicator for what the hold down criteria was (e.g., 800 or 1600)
T10_Logical_and_Continuous.criteria_indicator_All = criteria_indicator_All;
T10_Logical_and_Continuous.criteria_indicator_All = categorical(T10_Logical_and_Continuous.criteria_indicator_All);
%cumulative rewards earned
T10_Logical_and_Continuous.total_reward_indicator_All = total_reward_indicator_All;
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

%remove Nans for export to R
% T10_nansgone = T10_Logical_and_Continuous;
% T10_nansgone =T10_nansgone(~any(ismissing(T10_nansgone),2),:);
% writetable(T10_nansgone,'33-3-all-time.csv','QuoteStrings',true) %save table as
% csv
 %% This section is only applicable for optogenetic Stimulation, and is used to find variables relative to stimulation
 %e.g., press duration during/after stimulation to look for direct effect
 
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
% THIS REQUIRES THE quickENT FUNCTION FROM Timme and Lapish, 2018
 %duration entropy after stim/nostim
 lp_stim_n1_ent =quickEnt(lp_stim_n1);
 lp_no_stim_n1_ent = quickEnt(lp_no_stim_n1);
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
mouse_stim_n1_ent = [];
mouse_no_stim_n1_ent = [];
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
  % THIS REQUIRES THE quickENT FUNCTION FROM Timme and Lapish, 2018
    mouse_stim_n1_ent = [mouse_stim_n1_ent quickEnt(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_stim_idx_n1))];
    mouse_no_stim_n1_ent = [mouse_no_stim_n1_ent quickEnt(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_no_stim_idx_n1))];
    
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
mouse_stim_n1_ent = reshape(mouse_stim_n1_ent,[total_mice,day_number]);
mouse_no_stim_n1_ent = reshape(mouse_no_stim_n1_ent,[total_mice,day_number]);
mouse_means = reshape(mouse_means,[total_mice,day_number]);
mouse_std = reshape(mouse_std,[total_mice,day_number]);
mouse_mean_over_std =reshape(mouse_mean_over_std,[total_mice,day_number]);
mouse_med=reshape(mouse_med,[total_mice,day_number]);
mouse_iqr=reshape(mouse_iqr,[total_mice,day_number]);
mouse_iqr_over_med=reshape(mouse_iqr_over_med,[total_mice,day_number]);

%is there a correlation between predictability and task performance?
% figure()
% plot(correctCIprop_criteria_All(:,1)*100,correctCIprop_criteria_All(:,2))
% corr(correctCIprop_criteria_All(:,1)*100,correctCIprop_criteria_All(:,2));
table_of_preds = array2table(correctCIprop_criteria_All,'VariableNames',{'Correct_Pred','Criteria'});
table_of_preds.Correct_Pred = table_of_preds.Correct_Pred*100;
% fitlm(table_of_preds, 'Correct_Pred ~ Criteria');


%individual Regression predictions and criteria% by mouse x day
correctCIprop_All_by_day = reshape(correctCIprop_criteria_All(:,1)*100,[total_mice,day_number]);
crit_percent_All_by_day = reshape(correctCIprop_criteria_All(:,2),[total_mice,day_number]);
%get the means for each mouse (ie collapse days)
mouse_avg_predict = mean(correctCIprop_All_by_day,2);
mouse_avg_criteria = nanmean(crit_percent_All_by_day,2);
%plot the prediction vs criteria for mouse avg.
% plot(mouse_avg_criteria,mouse_avg_predict);
% corr(mouse_avg_criteria,mouse_avg_predict);

%complex model any better, or more correlated?
% figure()
% plot(correctCIprop_criteria_complex_All(:,1)*100,correctCIprop_criteria_complex_All(:,2));
% corr(correctCIprop_criteria_complex_All(:,1)*100,correctCIprop_criteria_complex_All(:,2));
% figure()
% corrplot(table_of_preds);

table_of_preds_complex = array2table(correctCIprop_criteria_complex_All,'VariableNames',{'Correct_Pred','Criteria'});
table_of_preds_complex.Correct_Pred = table_of_preds_complex.Correct_Pred*100;
% fitlm(table_of_preds_complex, 'Correct_Pred ~ Criteria');

%now use r^2 instead of prediction accuracy
% figure()
% plot(r_squared_criteria_All(:,1),r_squared_criteria_All(:,2))
table_of_r2 =  array2table(r_squared_criteria_All,'VariableNames',{'R2','Criteria'});
% fitlm(table_of_r2, 'R2 ~ Criteria');

r2_All_by_day = reshape(r_squared_criteria_All(:,1),[total_mice,day_number]);
% corr(mean(r2_All_by_day,2),mouse_avg_criteria/100);
% figure()
% plot(mean(r2_All_by_day,2),mouse_avg_criteria/100)

%  figure()
%  plot(r_squared_complex_criteria_All(:,1),r_squared_complex_criteria_All(:,2))
table_of_r2_complex =  array2table(r_squared_complex_criteria_All,'VariableNames',{'R2','Criteria'});
% fitlm(table_of_r2_complex, 'R2 ~ Criteria');


%is there a correlation between n-1 coef value and r2?
% figure()
% plot(indiv_lme_10_indiv_coef_All(2,:),r_squared_complex_criteria_All(:,2))
table_of_r2.Coef = indiv_lme_10_indiv_coef_All(2,:)';
% fitlm(table_of_r2, 'Coef~ Criteria ');

table_of_r2.Predictions = correctCIprop_criteria_All(:,1);

%add in mouse number
table_of_r2.Mouse = mouse_number_for_corrs;
%save to a csv for export to R for use in RMCORR
% writetable(table_of_r2,'m2_dms_sham_diffpreds_simple.csv','QuoteStrings',true) %save table as

% figure()
% corrplot(table_of_r2,'type','Pearson','testR','on')

table_of_r2_complex.Pred = correctCIprop_criteria_complex_All(:,1);


%add in mouse number
table_of_r2_complex.Mouse = mouse_number_for_corrs;
table_of_r2_complex.Day = day_number_for_corrs;

%save to a csv for export to R for use in RMCORR
% writetable(table_of_r2_complex,'m2_dms_sham_diffpreds_complex.csv','QuoteStrings',true) %save table as

table_of_r2_complex.time =total_time_and_rewards_All(:,1);
table_of_r2_complex.rew =total_time_and_rewards_All(:,2);
% table_of_r2_complex.Mouse =categorical(table_of_r2_complex.Mouse);
% table_of_r2_complex.Day =categorical(table_of_r2_complex.Day);

%lets create an lme where we can allow intercept and slope to vary across
% %mice
% 'R2 ~ Criteria + (1+Criteria|Mouse)';
% lme_r2_crit_complex_slopeint = fitlme(table_of_r2_complex,'R2 ~ Criteria*Day + (Day|Mouse)');
% plotPartialDependence(lme_r2_crit_complex_slopeint,{'Criteria','Mouse'})
% [B Bnames Bstats] =randomEffects(lme_r2_crit_complex_slopeint)
% [~,~,stats] = covarianceParameters(lme_r2_crit_complex_slopeint)
% plot(B(1:2:end),B(2:2:end),'r*')
% lme_r2_crit_complex_int_only = fitlme(table_of_r2_complex,'R2 ~ Criteria + (1|Mouse)');
% 
% compare(lme_r2_crit_complex_int_only,lme_r2_crit_complex_slopeint,'CheckNesting',true);

% figure()
% table_of_r2_complex.rewme = indiv_lme_reduced_coef_All(27,:)';
% 
% corrplot(table_of_r2_complex,'type','Pearson','testR','on')

 %now we will run correlations between Criteria% and CPD for n1 in the
 %simple and complex models, and the mov avg in the complex model
 table_of_CPD=  array2table(CPD_n1_criteria_simple_All,'VariableNames',{'N1_simple','Criteria'});
table_of_CPD.N1_complex = CPD_n1_criteria_complex_All(:,1);
table_of_CPD.MA_complex = CPD_MA_criteria_complex_All(:,1);
table_of_CPD.Mouse = mouse_number_for_corrs;
table_of_CPD.reward = CPD_n1rew_criteria_complex_All(:,1);

% table_of_CPD.rewme = indiv_lme_reduced_coef_All(14,:)';
table_of_CPD. r2 =r_squared_complex_criteria_All(:,1);
% 
% table_of_CPD.time =total_time_and_rewards_All(:,1);
% table_of_CPD.rew =total_time_and_rewards_All(:,2);
% 
% figure()
corrplot(table_of_CPD,'type','Pearson','testR','on')

% writetable(table_of_CPD,'M2_DMS_Lesion.csv','QuoteStrings',true) %save table as

%return the indicator variables to categorical for use in LMEs
T10_Logical_and_Continuous.day_indicator = categorical(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =categorical(T10_Logical_and_Continuous.mouse_indicator);
T10_Logical_and_Continuous.stim_indicator_All = categorical(T10_Logical_and_Continuous.stim_indicator_All);
T10_Logical_and_Continuous.stim_indicator_n1_All = categorical(T10_Logical_and_Continuous.stim_indicator_n1_All);
T10_Logical_and_Continuous.stim_indicator_n2_All = categorical(T10_Logical_and_Continuous.stim_indicator_n2_All);
T10_Logical_and_Continuous.stim_indicator_n3_All = categorical(T10_Logical_and_Continuous.stim_indicator_n3_All);
T10_Logical_and_Continuous.stim_indicator_n4_All = categorical(T10_Logical_and_Continuous.stim_indicator_n4_All);
T10_Logical_and_Continuous.stim_indicator_n5_All = categorical(T10_Logical_and_Continuous.stim_indicator_n5_All);
T10_Logical_and_Continuous.stim_indicator_n6_All = categorical(T10_Logical_and_Continuous.stim_indicator_n6_All);

 %% calculate entropy of all lps
 % THIS REQUIRES THE quickENT FUNCTION FROM Timme and Lapish, 2018
 lpent = quickEnt(T10_Logical_and_Continuous.LP_Durations_All);

%% shuffled order regressions
%now we need to make LMEs where we compare the shuffled with the real
%variable for each variable
%first the cannonical LME, with all variables and relevant interactions
model_spec_FULL_n6 = 'LP_Durations_All ~ moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator + moving_average_lp_length_n7andback_All:total_reward_indicator_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All + moving_average_lp_length_n7andback_All:up_state_idx_n1_All  + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator + n_minus_one_Durations_All:total_reward_indicator_All +  n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All + n_minus_one_Durations_All:up_state_idx_n1_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All + up_state_idx_n1_All + n_minus_two_Durations_All:LP_Timestamps_All + n_minus_two_Durations_All:criteria_percent_indicator + n_minus_two_Durations_All:total_reward_indicator_All + n_minus_two_Durations_All:HE_n_2_Indicator_All +  ipi2:n_minus_two_Durations_All + n_minus_two_Durations_All:n_minus_two_All + n_minus_two_Durations_All:up_state_idx_n2_All + HE_n_2_Indicator_All + ipi2 + n_minus_two_All + up_state_idx_n2_All + n_minus_three_Durations_All:LP_Timestamps_All + n_minus_three_Durations_All:criteria_percent_indicator + n_minus_three_Durations_All:total_reward_indicator_All + n_minus_three_Durations_All:HE_n_3_Indicator_All +  ipi3:n_minus_three_Durations_All + n_minus_three_Durations_All:n_minus_three_All + n_minus_three_Durations_All:up_state_idx_n3_All + HE_n_3_Indicator_All + ipi3 + n_minus_three_All + up_state_idx_n3_All + n_minus_four_Durations_All:LP_Timestamps_All + n_minus_four_Durations_All:criteria_percent_indicator + n_minus_four_Durations_All:total_reward_indicator_All + n_minus_four_Durations_All:HE_n_4_Indicator_All +  ipi4:n_minus_four_Durations_All + n_minus_four_Durations_All:n_minus_four_All + n_minus_four_Durations_All:up_state_idx_n4_All  + HE_n_4_Indicator_All + ipi4 + n_minus_four_All + up_state_idx_n4_All +n_minus_five_Durations_All:LP_Timestamps_All + n_minus_five_Durations_All:criteria_percent_indicator + n_minus_five_Durations_All:total_reward_indicator_All + n_minus_five_Durations_All:HE_n_5_Indicator_All +  ipi5:n_minus_five_Durations_All + n_minus_five_Durations_All:n_minus_five_All + n_minus_five_Durations_All:up_state_idx_n5_All  +HE_n_5_Indicator_All + ipi5 + n_minus_five_All + up_state_idx_n5_All +n_minus_six_Durations_All:LP_Timestamps_All + n_minus_six_Durations_All:criteria_percent_indicator + n_minus_six_Durations_All:total_reward_indicator_All + n_minus_six_Durations_All:HE_n_6_Indicator_All +  ipi6:n_minus_six_Durations_All + n_minus_six_Durations_All:n_minus_six_All + n_minus_six_Durations_All:up_state_idx_n6_All + HE_n_6_Indicator_All + ipi6 + n_minus_six_All + up_state_idx_n6_All + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + total_reward_indicator_All + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
lme_FULL_n6 = fitlme(T10_Logical_and_Continuous,model_spec_FULL_n6);
lme_FULL_n6_names = lme_FULL_n6.CoefficientNames';
lme_FULL_n6_coef = lme_FULL_n6.Coefficients.Estimate;
lme_FULL_n6_se = lme_FULL_n6.Coefficients.SE;
lme_FULL_n6_pval = lme_FULL_n6.Coefficients.pValue;
lme_FULL_n6_ftest = anova(lme_FULL_n6);
lme_FULL_n6_ftest_pval = lme_FULL_n6_ftest.pValue;
lme_FULL_n6_ftest_fstat = lme_FULL_n6_ftest.FStat	;

%now remove variables that were not significant in the max version
%BIC selected model
model_spec_reduced = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

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

%what about effects coding for categoricals?
lme_reduced_effects = fitlme(T10_Logical_and_Continuous,model_spec_reduced,'DummyVarCoding','effects');
lme_reduced_effects_coef = lme_reduced_effects.Coefficients.Estimate;
lme_reduced_effects_pval =  lme_reduced_effects.Coefficients.pValue;
anova(lme_reduced_effects);

%play with mouse interactions with n1 duration main effect 
model_spec_reduced_mouserew = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator) + (1|day_indicator) + (n_minus_one_Durations_All-1|mouse_indicator) ';

lme_reduced_mouserew = fitlme(T10_Logical_and_Continuous,model_spec_reduced_mouserew);
[b bnames bstats] =randomEffects(lme_reduced_mouserew);
% [r p] =corr(b(27:end,1), mouse_avg_criteria)
% plot(b(27:end,1), mouse_avg_criteria)
% corr(b(30:44,1), success_failure_mean_diff_overallmouse)
% plot(b(30:44,1), success_failure_mean_diff_overallmouse)
% corr(b(30:44,1), b(45:end,1))
% plot(b(30:44,1), b(45:end,1))
%find presses over 6000ms to potentially exlude
outl = find(lme_reduced.Residuals.Raw > 6000); 
lme_reduced_out = fitlme(T10_Logical_and_Continuous,model_spec_reduced,'Exclude',outl);
longdu = sort(T10_Logical_and_Continuous.LP_Durations_All);
%find presses over 10s and shorter than or equal to 20ms to potentialy
%exdlude
out2 =find(T10_Logical_and_Continuous.LP_Durations_All >= 10000 | T10_Logical_and_Continuous.LP_Durations_All <= 20);
lme_reduced_longshort = fitlme(T10_Logical_and_Continuous,model_spec_reduced,'Exclude',out2);
coef_norm_and_Exclude= [ lme_reduced.Coefficients.pValue	 lme_reduced_longshort.Coefficients.pValue];

%create the same for an lm for easy graphing
model_spec_reduce_lm = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator';
lm_reduced = fitlm(T10_Logical_and_Continuous,model_spec_reduce_lm);
% plotEffects(lm_reduced)
% figure()
% plotInteraction(lm_reduced,'ipi1','n_minus_one_Durations_All','predictions')
% figure()
% plotInteraction(lm_reduced,'ipi1','moving_average_lp_length_n7andback_All','predictions')
% lm_reduced_con=coefCI(lm_reduced);

%save simple model to fit to other datasets
% save('33_3lm.mat' ,'lm_reduced')


%create n10 back durations model continaing contorl vars for shuffling
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

model_spec_10_ma_totrew_mouse = 'LP_Durations_All~ criteria_percent_indicator+LP_Timestamps_All+ n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All + (1|mouse_indicator)+(1|day_indicator) + (n_minus_one_Durations_All-1|mouse_indicator)';
lme_10_ma_totrew_mouse = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew_mouse)
lme_10_ma_totrew_mouse_se = lme_10_ma_totrew_mouse.Coefficients.SE;
lme_10_ma_totrew_mouse_coef =  lme_10_ma_totrew_mouse.Coefficients.Estimate;
[b bnames bstats] =randomEffects(lme_10_ma_totrew_mouse);
% plot(b(1:12,1), b(27:end,1))
% corr(b(1:12,1), b(27:end,1))
% %is there a correlation between a mouse's avg percent correct and the
% %n-1 duration modifer? yes.
% plot(b(27:end,1), mouse_avg_criteria)
% corr(b(27:end,1), mouse_avg_criteria)



% temp_table = [b(1:15,1) b(45:end,1)];
% temp_table = array2table(temp_table,'VariableNames',{'MA','n1'})

% fitlm(temp_table,'n1~MA')
%does adding n10 improve an n9 model? 
model_spec_9_ma_totrew = 'LP_Durations_All~criteria_percent_indicator+LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
lme_9_ma_totrew = fitlme(T10_Logical_and_Continuous,model_spec_9_ma_totrew);
compare(lme_9_ma_totrew , lme_10_ma_totrew);
%Yes it does

%create n6 back durations model with control vars and MA
model_spec_6_ma_totrew = 'LP_Durations_All~ criteria_percent_indicator+LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+moving_average_lp_length_n7andback_All+(1|mouse_indicator)+(1|day_indicator)';
lme_6_ma_totrew = fitlme(T10_Logical_and_Continuous,model_spec_6_ma_totrew);
lme_6_ma_totrew_se = lme_6_ma_totrew.Coefficients.SE;
lme_6_ma_totrew_coef =  lme_6_ma_totrew.Coefficients.Estimate;

%now, within this reduced model we want to one at a time shuffle particular
%variables 
%this gets us all the data lumped together, with lps in order of occurrence
%(and order of mice/day like other variables) in rows and shuffles in
%columns
%the shuffled n-back for durations and logical are in a cell array that is
%number of shuffle rows by number of mouse/day columns
%each cell is number of n-back rows by number of total lps columns
% we need to make data tables that have, as columns, LP, n-back lp, day and
% mouse indicator. We then need to run a linear regression on this, store
% the coefficients, and repeat the process with the next shuffle

%% shuffling
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

%stimulation shuffle variables
shuf_stim_n_coef_all=[];
shuf_stim_n_SE_all=[];
shuf_stim_n1_coef_all = [];
shuf_stim_n1_SE_all = [];
shuf_stim_n2_coef_all = [];
shuf_stim_n2_SE_all = [];
shuf_stim_n3_coef_all = [];
shuf_stim_n3_SE_all = [];
shuf_stim_n4_coef_all = [];
shuf_stim_n4_SE_all = [];
shuf_stim_n5_coef_all = [];
shuf_stim_n5_SE_all = [];
shuf_stim_n6_coef_all = [];
shuf_stim_n6_SE_all = [];

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
%this is the reduced only, significant stuff model
% model_spec_reduced = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
% 
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

%% shuffled stimulus stuff
% 
% shuf_spec_stim_n = 'LP_Durations_All~shuf_stim_indicator_All*n_minus_one_Durations_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*shuf_stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n = fitlme(T_Both,shuf_spec_stim_n);
% shuf_stim_n_coef = shuf_stim_n.Coefficients.Estimate;
% shuf_stim_n_SE = shuf_stim_n.Coefficients.SE;
% shuf_stim_n_coef_all = [shuf_stim_n_coef_all shuf_stim_n_coef];
% shuf_stim_n_SE_all = [shuf_stim_n_SE_all shuf_stim_n_SE];
% 
% shuf_spec_stim_n1 = 'LP_Durations_All~n_minus_one_Durations_All*shuf_stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*shuf_stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n1 = fitlme(T_Both,shuf_spec_stim_n1);
% shuf_stim_n1_coef = shuf_stim_n1.Coefficients.Estimate;
% shuf_stim_n1_SE = shuf_stim_n1.Coefficients.SE;
% shuf_stim_n1_coef_all = [shuf_stim_n1_coef_all shuf_stim_n1_coef];
% shuf_stim_n1_SE_all = [shuf_stim_n1_SE_all shuf_stim_n1_SE];
% 
% shuf_spec_stim_n2 = 'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*shuf_stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n2 = fitlme(T_Both,shuf_spec_stim_n2);
% shuf_stim_n2_coef = shuf_stim_n2.Coefficients.Estimate;
% shuf_stim_n2_SE = shuf_stim_n2.Coefficients.SE;
% shuf_stim_n2_coef_all = [shuf_stim_n2_coef_all shuf_stim_n2_coef];
% shuf_stim_n2_SE_all = [shuf_stim_n2_SE_all shuf_stim_n2_SE];
% 
% shuf_spec_stim_n3 = 'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*shuf_stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n3 = fitlme(T_Both,shuf_spec_stim_n3);
% shuf_stim_n3_coef = shuf_stim_n3.Coefficients.Estimate;
% shuf_stim_n3_SE = shuf_stim_n3.Coefficients.SE;
% shuf_stim_n3_coef_all = [shuf_stim_n3_coef_all shuf_stim_n3_coef];
% shuf_stim_n3_SE_all = [shuf_stim_n3_SE_all shuf_stim_n3_SE];
% 
% shuf_spec_stim_n4 = 'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*shuf_stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n4 = fitlme(T_Both,shuf_spec_stim_n4);
% shuf_stim_n4_coef = shuf_stim_n4.Coefficients.Estimate;
% shuf_stim_n4_SE = shuf_stim_n4.Coefficients.SE;
% shuf_stim_n4_coef_all = [shuf_stim_n4_coef_all shuf_stim_n4_coef];
% shuf_stim_n4_SE_all = [shuf_stim_n4_SE_all shuf_stim_n4_SE];
% 
% shuf_spec_stim_n5 = 'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*shuf_stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n5 = fitlme(T_Both,shuf_spec_stim_n5);
% shuf_stim_n5_coef = shuf_stim_n5.Coefficients.Estimate;
% shuf_stim_n5_SE = shuf_stim_n5.Coefficients.SE;
% shuf_stim_n5_coef_all = [shuf_stim_n5_coef_all shuf_stim_n5_coef];
% shuf_stim_n5_SE_all = [shuf_stim_n5_SE_all shuf_stim_n5_SE];
% 
% shuf_spec_stim_n6 = 'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*shuf_stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)';
% shuf_stim_n6 = fitlme(T_Both,shuf_spec_stim_n6);
% shuf_stim_n6_coef = shuf_stim_n6.Coefficients.Estimate;
% shuf_stim_n6_SE = shuf_stim_n6.Coefficients.SE;
% shuf_stim_n6_coef_all = [shuf_stim_n6_coef_all shuf_stim_n6_coef];
% shuf_stim_n6_SE_all = [shuf_stim_n6_SE_all shuf_stim_n6_SE];

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
n_p_full_n_perm_all= NaN(length(lme_reduced_coef),1);

%duration me shuffles
n_p_full_n1_dur = sum( abs(shuf_full_n1dur_coef_all(2,:)) > abs(lme_reduced_coef(3)))/size(shuf_full_n1dur_coef_all,2);
n_p_full_n_perm_all(3) = n_p_full_n1_dur;

n_p_full_n2_dur = sum( abs(shuf_full_n2dur_coef_all(2,:)) > abs(lme_reduced_coef(4)))/size(shuf_full_n2dur_coef_all,2);
n_p_full_n_perm_all(4) = n_p_full_n2_dur;

n_p_full_n3_dur = sum( abs(shuf_full_n3dur_coef_all(2,:)) > abs(lme_reduced_coef(5)))/size(shuf_full_n3dur_coef_all,2);
n_p_full_n_perm_all(5) = n_p_full_n3_dur;

n_p_full_n4_dur = sum( abs(shuf_full_n4dur_coef_all(2,:)) > abs(lme_reduced_coef(6)))/size(shuf_full_n4dur_coef_all,2);
n_p_full_n_perm_all(6) = n_p_full_n4_dur;

n_p_full_n5_dur = sum( abs(shuf_full_n5dur_coef_all(2,:)) > abs(lme_reduced_coef(7)))/size(shuf_full_n5dur_coef_all,2);
n_p_full_n_perm_all(7) = n_p_full_n5_dur;

n_p_full_n6_dur = sum( abs(shuf_full_n6dur_coef_all(2,:)) > abs(lme_reduced_coef(8)))/size(shuf_full_n6dur_coef_all,2);
n_p_full_n_perm_all(8) = n_p_full_n6_dur;

n_p_full_ma_me = sum( abs(shuf_full_ma_coef_all(2,:)) > abs(lme_reduced_coef(14)))/size(shuf_full_ma_coef_all,2);
n_p_full_n_perm_all(14) = n_p_full_ma_me;
%also get mean and sem for the shuffled moving average for graphical purposes
full_ma_me_mean = mean(shuf_full_ma_coef_all(2,:));
full_ma_me_sem = std(shuf_full_ma_coef_all(2,:))/(sqrt(size( shuf_full_ma_coef_all(2,:),2)));

%now compare this to when we shuffle main effects and interactions
%individually
n_p_full_n_perm_all_indivshuffles =n_p_full_n_perm_all;

n_p_full_ipi_me = sum( abs(shuf_full_ipi_me_coef_all(2,:)) > abs(lme_reduced_coef(12)))/size(shuf_full_ipi_me_coef_all,2);
n_p_full_ipi_int = sum( abs(shuf_full_ipi_n1_coef_all(15,:)) > abs(lme_reduced_coef(19)))/size(shuf_full_ipi_n1_coef_all,2);
n_p_full_ipi_MA = sum( abs(shuf_full_ipi_avg_coef_all(21,:)) > abs(lme_reduced_coef(25)))/size(shuf_full_ipi_avg_coef_all,2);
n_p_full_n_perm_all_indivshuffles(12) = n_p_full_ipi_me;
n_p_full_n_perm_all_indivshuffles(19) = n_p_full_ipi_int;
n_p_full_n_perm_all_indivshuffles(25) = n_p_full_ipi_MA;

n_p_full_ipi2_me = sum( abs(shuf_full_ipi2_me_coef_all(2,:)) > abs(lme_reduced_coef(13)))/size(shuf_full_ipi2_me_coef_all,2);
n_p_full_ipi2_int = sum( abs(shuf_full_ipi2_n2_coef_all(16,:)) > abs(lme_reduced_coef(20)))/size(shuf_full_ipi2_n2_coef_all,2);
n_p_full_n_perm_all_indivshuffles(13) = n_p_full_ipi2_me;
n_p_full_n_perm_all_indivshuffles(20) = n_p_full_ipi2_int;

n_p_full_he_me = sum( abs(shuf_full_he_me_coef_all(2,:)) > abs(lme_reduced_coef(10)))/size(shuf_full_he_me_coef_all,2);
n_p_full_he_int = sum( abs(shuf_full_he_n1_coef_all(15,:)) > abs(lme_reduced_coef(17)))/size(shuf_full_he_n1_coef_all,2);
n_p_full_he_MA = sum( abs(shuf_full_he_avg_coef_all(21,:)) > abs(lme_reduced_coef(23)))/size(shuf_full_he_avg_coef_all,2);
n_p_full_n_perm_all_indivshuffles(10) = n_p_full_he_me;
n_p_full_n_perm_all_indivshuffles(17) = n_p_full_he_int;
n_p_full_n_perm_all_indivshuffles(23) = n_p_full_he_MA;

n_p_full_rew_me = sum( abs(shuf_full_rew_me_coef_all(2,:)) > abs(lme_reduced_coef(2)))/size(shuf_full_rew_me_coef_all,2);
n_p_full_rew_int = sum( abs(shuf_full_rew_n1_coef_all(15,:)) > abs(lme_reduced_coef(15)))/size(shuf_full_rew_n1_coef_all,2);
n_p_full_rew_MA = sum( abs(shuf_full_rew_avg_coef_all(21,:)) > abs(lme_reduced_coef(21)))/size(shuf_full_rew_avg_coef_all,2);
n_p_full_n_perm_all_indivshuffles(2) = n_p_full_rew_me;
n_p_full_n_perm_all_indivshuffles(15) = n_p_full_rew_int;
n_p_full_n_perm_all_indivshuffles(21) = n_p_full_rew_MA;

n_p_full_ts_me = sum( abs(shuf_full_ts_me_coef_all(2,:)) > abs(lme_reduced_coef(11)))/size(shuf_full_ts_me_coef_all,2);
n_p_full_ts_int = sum( abs(shuf_full_ts_n1_coef_all(15,:)) > abs(lme_reduced_coef(18)))/size(shuf_full_ts_n1_coef_all,2);
n_p_full_ts_MA = sum( abs(shuf_full_ts_avg_coef_all(21,:)) > abs(lme_reduced_coef(24)))/size(shuf_full_ts_avg_coef_all,2);
n_p_full_n_perm_all_indivshuffles(11) = n_p_full_ts_me;
n_p_full_n_perm_all_indivshuffles(18) = n_p_full_ts_int;
n_p_full_n_perm_all_indivshuffles(24) = n_p_full_ts_MA;

n_p_full_crit_me = sum( abs(shuf_full_crit_me_coef_all(2,:)) > abs(lme_reduced_coef(9)))/size(shuf_full_crit_me_coef_all,2);
n_p_full_crit_int = sum( abs(shuf_full_crit_n1_coef_all(15,:)) > abs(lme_reduced_coef(16)))/size(shuf_full_crit_n1_coef_all,2);
n_p_full_crit_MA = sum( abs(shuf_full_crit_avg_coef_all(21,:)) > abs(lme_reduced_coef(22)))/size(shuf_full_crit_avg_coef_all,2);
n_p_full_n_perm_all_indivshuffles(9) = n_p_full_crit_me;
n_p_full_n_perm_all_indivshuffles(16) = n_p_full_crit_int;
n_p_full_n_perm_all_indivshuffles(22) = n_p_full_crit_MA;

%do it based on aic instead of coef value
n_p_full_n_perm_all_AIC= NaN(length(lme_reduced_coef),1);
lme_reduced_AIC = lme_reduced.ModelCriterion.AIC;

%with AIC, if the non-shuffled version is smaller, it is better. Thus we
%want to see how many times the shuffled version of AIC is larger than the
%acutal
%if AIC in the shuffled is smaller - that means it's "better" 
%that happens 12/100 times with n1. Thus, 12 of the shuffled ones we get a
%better fit to teh data, and should be rejecteed
n_p_full_n1_dur = sum( shuf_full_n1dur_AIC_all < lme_reduced_AIC)/size(shuf_full_n1dur_AIC_all,2);
n_p_full_n_perm_all_AIC(3) = n_p_full_n1_dur;

n_p_full_n2_dur = sum( abs(shuf_full_n2dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n2dur_AIC_all,2);
n_p_full_n_perm_all_AIC(4) = n_p_full_n2_dur;

n_p_full_n3_dur = sum( abs(shuf_full_n3dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n3dur_AIC_all,2);
n_p_full_n_perm_all_AIC(5) = n_p_full_n3_dur;

n_p_full_n4_dur = sum( abs(shuf_full_n4dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n4dur_AIC_all,2);
n_p_full_n_perm_all_AIC(6) = n_p_full_n4_dur;

n_p_full_n5_dur = sum( abs(shuf_full_n5dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n5dur_AIC_all,2);
n_p_full_n_perm_all_AIC(7) = n_p_full_n5_dur;

n_p_full_n6_dur = sum( abs(shuf_full_n6dur_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_n6dur_AIC_all,2);
n_p_full_n_perm_all_AIC(8) = n_p_full_n6_dur;

n_p_full_ma_me = sum( abs(shuf_full_ma_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ma_AIC_all,2);
n_p_full_n_perm_all_AIC(14) = n_p_full_ma_me;

n_p_full_ipi_me = sum( abs(shuf_full_ipi_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_me_AIC_all,2);
n_p_full_ipi_int = sum( abs(shuf_full_ipi_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_n1_AIC_all,2);
n_p_full_ipi_MA = sum( abs(shuf_full_ipi_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi_avg_AIC_all,2);
n_p_full_n_perm_all_AIC(12) = n_p_full_ipi_me;
n_p_full_n_perm_all_AIC(19) = n_p_full_ipi_int;
n_p_full_n_perm_all_AIC(25) = n_p_full_ipi_MA;

n_p_full_ipi2_me = sum( abs(shuf_full_ipi2_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi2_me_AIC_all,2);
n_p_full_ipi2_int = sum( abs(shuf_full_ipi2_n2_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ipi2_n2_AIC_all,2);
n_p_full_n_perm_all_AIC(13) = n_p_full_ipi2_me;
n_p_full_n_perm_all_AIC(20) = n_p_full_ipi2_int;

n_p_full_he_me = sum( abs(shuf_full_he_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_me_AIC_all,2);
n_p_full_he_int = sum( abs(shuf_full_he_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_n1_AIC_all,2);
n_p_full_he_MA = sum( abs(shuf_full_he_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_he_avg_AIC_all,2);
n_p_full_n_perm_all_AIC(10) = n_p_full_he_me;
n_p_full_n_perm_all_AIC(17) = n_p_full_he_int;
n_p_full_n_perm_all_AIC(23) = n_p_full_he_MA;

n_p_full_rew_me = sum( abs(shuf_full_rew_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_me_AIC_all,2);
n_p_full_rew_int = sum( abs(shuf_full_rew_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_n1_AIC_all,2);
n_p_full_rew_MA = sum( abs(shuf_full_rew_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_rew_avg_AIC_all,2);
n_p_full_n_perm_all_AIC(2) = n_p_full_rew_me;
n_p_full_n_perm_all_AIC(15) = n_p_full_rew_int;
n_p_full_n_perm_all_AIC(21) = n_p_full_rew_MA;

n_p_full_ts_me = sum( abs(shuf_full_ts_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_me_AIC_all,2);
n_p_full_ts_int = sum( abs(shuf_full_ts_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_n1_AIC_all,2);
n_p_full_ts_MA = sum( abs(shuf_full_ts_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_ts_avg_AIC_all,2);
n_p_full_n_perm_all_AIC(11) = n_p_full_ts_me;
n_p_full_n_perm_all_AIC(18) = n_p_full_ts_int;
n_p_full_n_perm_all_AIC(24) = n_p_full_ts_MA;

n_p_full_crit_me = sum( abs(shuf_full_crit_me_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_me_AIC_all,2);
n_p_full_crit_int = sum( abs(shuf_full_crit_n1_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_n1_AIC_all,2);
n_p_full_crit_MA = sum( abs(shuf_full_crit_avg_AIC_all) < abs(lme_reduced_AIC))/size(shuf_full_crit_avg_AIC_all,2);
n_p_full_n_perm_all_AIC(9) = n_p_full_crit_me;
n_p_full_n_perm_all_AIC(16) = n_p_full_crit_int;
n_p_full_n_perm_all_AIC(22) = n_p_full_crit_MA;
     

%% Predictions and diagnostics of the behavioral model
%it turns out that all the variables selected via AIC are also selected by
%the permutation test. Thus, we will keep the "model_spec_reduced" model
%here
model_spec_reduced = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

model_spec_reduced_nocrit = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All +  n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All +  (1|mouse_indicator)+(1|day_indicator)';
lme_reduced_nocrit = fitlme(T10_Logical_and_Continuous,model_spec_reduced_nocrit);
% lme_reduced
% lm_reduced
%Save the simple linear regression. Can be used to see how well other
%datasets can be predicted 
% save('simplelm.mat' ,'lm_reduced')

% plotEffects(lm_reduced)

%try to replicate the lm ploteffects graph with the lme
%take the main effect coef and ci at 0, then multiply both by a reasonable
%number (e.g., 1000ms for lp length), or just take the et and cof interval
%for categorical

% labels = lme_reduced.CoefficientNames;
% coefestimates = lme_reduced.Coefficients.Estimate;
% plotAdjustedResponse(lm_reduced,'n_minus_one_Durations_All')

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
 
 %% Individual Mouse Data for graphical purposes
 %Get some individual mice to show how well the model predicts individual
 %data/for graphing individual data
 %looking at the final day of training

% 40523 - 40809 mouse 8 (actual mouse 12)
% 40928 - 41058 mouse 10 (actual mouse 14 in 33-3, as mice are labeled in nummerical order for the matlab script)
% 41260 - end  mouse 12 (actual mouse 16)
mouse10_lps = T10_Logical_and_Continuous.LP_Durations_All(40928:41058);
mouse10rew =T10_Logical_and_Continuous.LP_All(40928:41058);
mouse10lps_smooth = smoothdata(mouse10_lps,'gaussian',5);

mouse8lps = T10_Logical_and_Continuous.LP_Durations_All(40523:40809);
mouse8lpspred = yhat((40523:40809 ));
mouse8rew =T10_Logical_and_Continuous.LP_All(40523:40809);
mouse8lps_smooth = smoothdata(mouse8lps,'gaussian',5);

mouse12lps_smooth = smoothdata(T10_Logical_and_Continuous.LP_Durations_All(41260:end),'gaussian',[10 0]);
mouse12lps_shuf_smooth = smoothdata(Shuffled_Durations_All(41260:end),'gaussian',[10 0]);
correctish_pred =T10_Logical_and_Continuous.LP_Durations_All(41260:end) <=yhatCI(41260:end,2) &  T10_Logical_and_Continuous.LP_Durations_All(41260:end) >=yhatCI(41260:end,1);
correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.LP_Durations_All(41260:end));

%% Can we predict specific mice or days better than others?
%use a model without criteria% in it. 
  [yhat yhatCI yhatDF]= predict(lme_reduced_nocrit,'Simultaneous',true);
%     [yhat yhatCI yhatDF]= predict(lme_reduced_nocrit,'Simultaneous',false);

  lme_reduced_nocrit_residuals = residuals(lme_reduced_nocrit);

it_d =length(unique(T10_Logical_and_Continuous.day_indicator));
it_m = length(unique(T10_Logical_and_Continuous.mouse_indicator));
T10_Logical_and_Continuous.day_indicator =double(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =double(T10_Logical_and_Continuous.mouse_indicator);
correctCI_prop_All=[];
Residuals_by_mouse =[];
rew_red_actual_by_mouse =[];
for day_mean = 1:it_d 
     day_idx = find(T10_Logical_and_Continuous.day_indicator == day_mean);
   
for mouse_mean = 1:it_m
    
    mouse_idx = find(T10_Logical_and_Continuous.mouse_indicator ==mouse_mean);
    mouse_by_day_idx = intersect(day_idx,mouse_idx);
    correctish_pred =T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx) <=yhatCI(mouse_by_day_idx,2) &  T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx) >=yhatCI(mouse_by_day_idx,1);
    correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.LP_Durations_All(mouse_by_day_idx));
    correctCI_prop_All = [correctCI_prop_All correctCI_prop];
    %look at residuals by mouse
    Residuals_by_mouse =[Residuals_by_mouse nanmean(lme_reduced_nocrit_residuals(mouse_by_day_idx))];
    %look at prediction of successfailure
    rew_red_actual_by_mouse = [rew_red_actual_by_mouse sum(rew_red_actual(mouse_by_day_idx,1) == rew_red_actual(mouse_by_day_idx,2))/length(mouse_by_day_idx) ];
    
    
end

end
x = rew_red_actual(mouse_by_day_idx,:)

%correlation between predictability and criteria%?
corr(correctCIprop_criteria_All(:,2));
table_overall_preds = array2table(correctCI_prop_All', 'VariableNames',{'Predictions'});
table_overall_preds.residuals = Residuals_by_mouse';
table_overall_preds.crit =correctCIprop_criteria_All(:,2);
table_overall_preds.reward = rew_red_actual_by_mouse';
table_overall_preds.mouse =mouse_number_for_corrs;
table_overall_preds.day =day_number_for_corrs;
corrplot(table_overall_preds,'testR','on');
corrplot(table_overall_preds(1:72,:),'testR','on');
corrplot(table_overall_preds(73:end,:),'testR','on');

%mouse x days of training proportion predicted correctly
correctCI_prop_All_by_session=reshape(correctCI_prop_All,[total_mice,day_number]);
%lets get both means by mouse and by day
correctCI_prop_All_mean_by_mouse = mean(correctCI_prop_All_by_session,2);
correctCI_prop_All_mean_by_day = mean(correctCI_prop_All_by_session,1);
crit_mean_by_mouse = mean(crit_percent_All_by_day,2);
correctCIprop_criteria_All(:,2);

corr(crit_mean_by_mouse,correctCI_prop_All_mean_by_mouse);
% writetable(table_overall_preds,'33_3_Overall_Preds_nocrit.csv','QuoteStrings',true)
% %save table for rmcorr in R

% overall_crit_pred_table = array2table(correctCI_prop_All_mean_by_mouse, 'VariableNames',{'Predictions'}];
% corr(correctCI_prop_All_mean_by_mouse,crit_mean_by_mouse)
%

T10_Logical_and_Continuous.day_indicator =categorical(T10_Logical_and_Continuous.day_indicator);
T10_Logical_and_Continuous.mouse_indicator =categorical(T10_Logical_and_Continuous.mouse_indicator);

%% Does adding "future" decisions improve predictions?
%create a model that predicts the middle given 5 beore and 5 after
model_spec_history_and_future= 'n_minus_five_Durations_All~LP_Durations_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All +criteria_percent_indicator+(1|mouse_indicator)+(1|day_indicator)';
lme_history_and_future = fitlme(T10_Logical_and_Continuous,model_spec_history_and_future);
 [yhat yhatCI yhatDF]= predict(lme_history_and_future,T10_Logical_and_Continuous,'Simultaneous',true);
   correctish_pred =T10_Logical_and_Continuous.n_minus_five_Durations_All <=yhatCI(:,2) &  T10_Logical_and_Continuous.n_minus_five_Durations_All >=yhatCI(:,1);
 correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.n_minus_five_Durations_All)  ;
   %15.2% for simple model with n5 prior and n5 post vs. the 18.1%
   %predicted accurately in the n10 prior model
   %suggests that mice are using history, and that it is not simply due to
   %relatedness of presses 


  
 x= 1:length(T10_Logical_and_Continuous.LP_Durations_All);
 x=x';
 y = T10_Logical_and_Continuous.LP_Durations_All;
 scatter(x,y)   % the data circles
hold on
plot(x,yhat,':')  % the regression line
arrayfun(@(x,y,yhat) line([x x],[y yhat],'linestyle','-'),x,y,yhat)  % the stem lines
 
smootheddata = smoothdata(T10_Logical_and_Continuous.LP_Durations_All,'sgolay',10);
smoothpred = smoothdata(yhat,'sgolay',10);
T10_Logical_and_Continuous.LP_All = double(T10_Logical_and_Continuous.LP_All);
rewindices = find(T10_Logical_and_Continuous.LP_All ==1);
T10_Logical_and_Continuous.LP_All = categorical(T10_Logical_and_Continuous.LP_All);



%% Playing around with Random effect specification
% Does adding a random effect for criteria improve the model?
model_spec_reduce_critint= 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)+ (1|criteria_indicator_All)';
lme_reduce_critint = fitlme(T10_Logical_and_Continuous,model_spec_reduce_critint);
[B,Bnames,stats]  = randomEffects(lme_reduce_critint);
lme_reduced.ModelCriterion.AIC; - lme_reduce_critint.ModelCriterion.AIC;

pt1 = linspace(min(T10_Logical_and_Continuous.moving_average_lp_length_n7andback_All),max(T10_Logical_and_Continuous.moving_average_lp_length_n7andback_All),50)';
pt2 = linspace(min(T10_Logical_and_Continuous.ipi1),max(T10_Logical_and_Continuous.ipi1),50)';
% 
% plotPartialDependence(lme_reduce_critint,{'moving_average_lp_length_n7andback_All','ipi1'},'QueryPoints',[pt1 pt2])
% plotPartialDependence(lme_reduce_critint,{'moving_average_lp_length_n7andback_All','n_minus_one_All'})

MdlStd = fitrsvm(T10_Logical_and_Continuous,'LP_Durations_All','KernelFunction','gaussian','KernelScale','auto','Standardize',true)

model_spec_reduce_critint_crit_nest= 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All + LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator) + (1|day_indicator:criteria_indicator_All)';
lme_reduce_critint_critnest = fitlme(T10_Logical_and_Continuous,model_spec_reduce_critint_crit_nest);
[B,Bnames,stats]  = randomEffects(lme_reduce_critint_critnest);

%% Stimulation LMEs

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

%add in he ints
lme3_time_stim_n_int_ctl_ma_int_full_hen = fitlme(T10_Logical_and_Continuous,'LP_Durations_All ~HE_n_1_Indicator_All:stim_indicator_All + HE_n_1_Indicator_All:stim_indicator_All:n_minus_one_Durations_All + stim_indicator_All*n_minus_one_Durations_All + stim_indicator_n1_All*n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All + criteria_percent_indicator+ipi2*n_minus_two_Durations_All+n_minus_one_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n_1_Indicator_All+LP_Timestamps_All*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All + ipi1:moving_average_lp_length_n7andback_All +  moving_average_lp_length_n7andback_All:n_minus_one_All+moving_average_lp_length_n7andback_All:criteria_percent_indicator  + LP_Timestamps_All:moving_average_lp_length_n7andback_All+stim_indicator_All*moving_average_lp_length_n7andback_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All  + (1|mouse_indicator)+(1|day_indicator)');
compare(lme3_time_stim_n_int_ctl_ma_int_full,lme3_time_stim_n_int_ctl_ma_int_full_hen)
lme3_time_stim_n_int_ctl_ma_int_full_hen1 = fitlme(T10_Logical_and_Continuous,'LP_Durations_All ~HE_n_1_Indicator_All:stim_indicator_n1_All + HE_n_1_Indicator_All:stim_indicator_n1_All:n_minus_one_Durations_All + stim_indicator_All*n_minus_one_Durations_All + stim_indicator_n1_All*n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All + criteria_percent_indicator+ipi2*n_minus_two_Durations_All+n_minus_one_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n_1_Indicator_All+LP_Timestamps_All*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+criteria_percent_indicator:n_minus_one_Durations_All+moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All + ipi1:moving_average_lp_length_n7andback_All +  moving_average_lp_length_n7andback_All:n_minus_one_All+moving_average_lp_length_n7andback_All:criteria_percent_indicator  + LP_Timestamps_All:moving_average_lp_length_n7andback_All+stim_indicator_All*moving_average_lp_length_n7andback_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All  + (1|mouse_indicator)+(1|day_indicator)');
compare(lme3_time_stim_n_int_ctl_ma_int_full,lme3_time_stim_n_int_ctl_ma_int_full_hen1)

%% stimulation lmes
% lme3_time_stim = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_coef = lme3_time_stim.Coefficients.Estimate;
% lme3_time_stim_SE = lme3_time_stim.Coefficients.SE;
% lme3_time_stim_name = lme3_time_stim.Coefficients.Name;
% lme3_time_stim_pval= lme3_time_stim.Coefficients.pValue;
% % 
% % lme3_time_stim_lm = fitlm(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All');
% % figure()
% % plotInteraction(lme3_time_stim_lm,'n_minus_one_Durations_All','stim_indicator_n1_All')
% 
% 
% lme3_time_stim_n = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~stim_indicator_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_n_int = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~stim_indicator_All*n_minus_one_Durations_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_n_int_name = lme3_time_stim_n_int.Coefficients.Name;
% lme3_time_stim_n_int_coef = lme3_time_stim_n_int.Coefficients.Estimate;
% lme3_time_stim_n_int_SE = lme3_time_stim_n_int.Coefficients.SE;
% lme3_time_stim_n_int_pval = lme3_time_stim_n_int.Coefficients.pValue;
% % lme3_time_stim_coef
% 
% %add in all the non lp modifiers (criteria%, total reward, MA)
% lme3_time_stim_n_int_ctl_ma = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~stim_indicator_All*n_minus_one_Durations_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)');
% %add in interaction with n and n1 with moving average
% lme3_time_stim_n_int_ctl_ma_int = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+moving_average_lp_length_n7andback_All*stim_indicator_All+moving_average_lp_length_n7andback_All*stim_indicator_n1_All+criteria_percent_indicator+total_reward_indicator_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_n_int_ctl_ma_int_coef = lme3_time_stim_n_int_ctl_ma_int.Coefficients.Estimate;
% 
% 
% lme3_time_stim_ipi = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All +stim_indicator_n1_All*ipi1+ n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_re = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*stim_indicator_n1_All +stim_indicator_n1_All*n_minus_one_All*n_minus_one_Durations_All+ n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_ipi2 = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All:ipi1+n_minus_one_Durations_All*stim_indicator_n1_All +stim_indicator_n1_All*ipi1+ n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_time_stim_ipi3 = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All:ipi1:stim_indicator_n1_All+n_minus_one_Durations_All:ipi1+n_minus_one_Durations_All*stim_indicator_n1_All +stim_indicator_n1_All*ipi1+ n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% 
% lme3_time_stim_ma = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~moving_mean*stim_indicator_n1_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% 
% lme3_time_stim_HE = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All*HE_n_1_Indicator_All+n_minus_one_Durations_All:stim_indicator_n1_All:HE_n_1_Indicator_All+n_minus_one_Durations_All*stim_indicator_n1_All +HE_n_1_Indicator_All+ n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% 
% lme3_time_stim_state = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~up_state_idx_n1_All:n_minus_one_Durations_All:stim_indicator_n1_All+up_state_idx_n1_All*n_minus_one_Durations_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All*stim_indicator_n2_All + n_minus_three_Durations_All*stim_indicator_n3_All + n_minus_four_Durations_All*stim_indicator_n4_All+n_minus_five_Durations_All*stim_indicator_n5_All+n_minus_six_Durations_All*stim_indicator_n6_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% 
% %model with all interactions up to n6 - add in n1 stim
% model_spec_n6_allints_ma_stim = 'LP_Durations_All ~n_minus_one_Durations_All*stim_indicator_n1_All+ n_minus_one_Durations_All*HE_n_1_Indicator_All + ipi1:n_minus_one_Durations_All +n_minus_one_Durations_All:up_state_idx_n1_All + n_minus_one_Durations_All:n_minus_one_All + n_minus_two_Durations_All + n_minus_two_Durations_All:HE_n_2_Indicator_All + ipi2:n_minus_two_Durations_All +n_minus_two_Durations_All:up_state_idx_n2_All + n_minus_two_Durations_All:n_minus_two_All + n_minus_three_Durations_All + n_minus_three_Durations_All:HE_n_3_Indicator_All + ipi3:n_minus_three_Durations_All +n_minus_three_Durations_All:up_state_idx_n3_All + n_minus_three_Durations_All:n_minus_three_All + n_minus_four_Durations_All + n_minus_four_Durations_All:HE_n_4_Indicator_All + ipi4:n_minus_four_Durations_All +n_minus_four_Durations_All:up_state_idx_n4_All + n_minus_four_Durations_All:n_minus_four_All + n_minus_five_Durations_All + n_minus_five_Durations_All:HE_n_5_Indicator_All + ipi5:n_minus_five_Durations_All +n_minus_five_Durations_All:up_state_idx_n5_All + n_minus_five_Durations_All:n_minus_five_All + n_minus_six_Durations_All + n_minus_six_Durations_All:HE_n_6_Indicator_All + ipi6:n_minus_six_Durations_All +n_minus_six_Durations_All:up_state_idx_n6_All + n_minus_six_Durations_All:n_minus_six_All+moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All: HE_n_1_Indicator_All + moving_average_lp_length_n7andback_All:ipi1 + moving_average_lp_length_n7andback_All:up_state_idx_n1_All + moving_average_lp_length_n7andback_All:n_minus_one_All + (1|mouse_indicator)+(1|day_indicator)';
% lme_n6_allints_ma_stim = fitlme(T10_Logical_and_Continuous,model_spec_n6_allints_ma_stim);
% 
% %compaire adding in stim to a model with only durations
% lme3_n6 = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% lme3_n6_s1me = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All+stim_indicator_n1_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% compare(lme3_n6,lme3_n6_s1me);
% lme3_n6_s1me_int = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All+n_minus_one_Durations_All*stim_indicator_n1_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% compare(lme3_n6_s1me,lme3_n6_s1me_int);
% 
% lme3_n6_s1me_int_nly = fitlme(T10_Logical_and_Continuous,'LP_Durations_All~n_minus_one_Durations_All+n_minus_one_Durations_All:stim_indicator_n1_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All +LP_Timestamps_All+(1|mouse_indicator)+(1|day_indicator)');
% compare(lme3_n6,lme3_n6_s1me_int_nly);

%% means and sem for stim
% shuf_stim_n_coef_mean = mean(shuf_stim_n_coef_all,2);
% shuf_stim_n_coef_std = std(shuf_stim_n_coef_all,0,2);
% shuf_stim_n_coef_sem = shuf_stim_n_coef_std/(sqrt(size(shuf_stim_n_coef_all,2)));
% 
% shuf_stim_n1_coef_mean = mean(shuf_stim_n1_coef_all,2);
% shuf_stim_n1_coef_std = std(shuf_stim_n1_coef_all,0,2);
% shuf_stim_n1_coef_sem = shuf_stim_n1_coef_std/(sqrt(size(shuf_stim_n1_coef_all,2)));
% 
% shuf_stim_n2_coef_mean = mean(shuf_stim_n2_coef_all,2);
% shuf_stim_n2_coef_std = std(shuf_stim_n2_coef_all,0,2);
% shuf_stim_n2_coef_sem = shuf_stim_n2_coef_std/(sqrt(size(shuf_stim_n2_coef_all,2)));
% 
% shuf_stim_n3_coef_mean = mean(shuf_stim_n3_coef_all,2);
% shuf_stim_n3_coef_std = std(shuf_stim_n3_coef_all,0,2);
% shuf_stim_n3_coef_sem = shuf_stim_n3_coef_std/(sqrt(size(shuf_stim_n3_coef_all,2)));
% 
% shuf_stim_n4_coef_mean = mean(shuf_stim_n4_coef_all,2);
% shuf_stim_n4_coef_std = std(shuf_stim_n4_coef_all,0,2);
% shuf_stim_n4_coef_sem = shuf_stim_n4_coef_std/(sqrt(size(shuf_stim_n4_coef_all,2)));
% 
% shuf_stim_n5_coef_mean = mean(shuf_stim_n5_coef_all,2);
% shuf_stim_n5_coef_std = std(shuf_stim_n5_coef_all,0,2);
% shuf_stim_n5_coef_sem = shuf_stim_n5_coef_std/(sqrt(size(shuf_stim_n5_coef_all,2)));
% 
% shuf_stim_n6_coef_mean = mean(shuf_stim_n6_coef_all,2);
% shuf_stim_n6_coef_std = std(shuf_stim_n6_coef_all,0,2);
% shuf_stim_n6_coef_sem = shuf_stim_n6_coef_std/(sqrt(size(shuf_stim_n6_coef_all,2)));

% % permutation tests for stim
% n_p_stim_n = sum( abs(shuf_stim_n_coef_all(18,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(12)))/size(shuf_stim_n_coef_all,2);
% n_p_stim_n_int = sum( abs(shuf_stim_n_coef_all(26,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(19)))/size(shuf_stim_n_coef_all,2);
% n_p_stim_n_MA = sum( abs(shuf_stim_n_coef_all(27,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(20)))/size(shuf_stim_n_coef_all,2);
% 
% n_p_stim_n1 = sum( abs(shuf_stim_n1_coef_all(18,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(13)))/size(shuf_stim_n1_coef_all,2);
% n_p_stim_n1_int = sum( abs(shuf_stim_n1_coef_all(26,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(21)))/size(shuf_stim_n1_coef_all,2);
% n_p_stim_n1_MA = sum( abs(shuf_stim_n1_coef_all(27,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(22)))/size(shuf_stim_n1_coef_all,2);
% 
% n_p_stim_n2 = sum( abs(shuf_stim_n2_coef_all(15,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(11)))/size(shuf_stim_n2_coef_all,2);
% n_p_stim_n2_int = sum( abs(shuf_stim_n2_coef_all(22,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(18)))/size(shuf_stim_n2_coef_all,2);
% 
% n_p_stim_n3 = sum( abs(shuf_stim_n3_coef_all(15,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(12)))/size(shuf_stim_n3_coef_all,2);
% n_p_stim_n3_int = sum( abs(shuf_stim_n3_coef_all(22,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(19)))/size(shuf_stim_n3_coef_all,2);
% 
% n_p_stim_n4 = sum( abs(shuf_stim_n4_coef_all(15,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(13)))/size(shuf_stim_n4_coef_all,2);
% n_p_stim_n4_int = sum( abs(shuf_stim_n4_coef_all(22,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(20)))/size(shuf_stim_n4_coef_all,2);
% 
% n_p_stim_n5 = sum( abs(shuf_stim_n5_coef_all(15,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(14)))/size(shuf_stim_n5_coef_all,2);
% n_p_stim_n5_int = sum( abs(shuf_stim_n5_coef_all(22,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(21)))/size(shuf_stim_n5_coef_all,2);
% 
% n_p_stim_n6 = sum( abs(shuf_stim_n6_coef_all(15,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(15)))/size(shuf_stim_n6_coef_all,2);
% n_p_stim_n6_int = sum( abs(shuf_stim_n6_coef_all(22,:)) > abs(lme3_time_stim_n_int_ctl_ma_int_coef(22)))/size(shuf_stim_n6_coef_all,2);

 %% for probability add in main effect and interaction of success + rew vs success no reward
% remove n_minus_one_all effects, since it only cares whether duration was long enough, not whether there was reward
% (this only matters if reward is probbalistic)
% instead use the reward_indicator_n1_All variable
% where 0 = fail, no reward, 1 = success, no reward, and 2 = success, yes
% reward
T10_Logical_and_Continuous.reward_indicator_n1_All = categorical(T10_Logical_and_Continuous.reward_indicator_n1_All);
T10_Logical_and_Continuous.reward_indicator_All = categorical(T10_Logical_and_Continuous.reward_indicator_All);

T10_Logical_and_Continuous.total_reward_indicator_prob_All = total_reward_indicator_prob_All;

%BIC selected model
model_spec_moderate_remove_proability = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

%BIC selected model + cumulative success/cumulative reward
model_spec_moderate_remove_proability = 'LP_Durations_All ~ total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All + total_reward_indicator_All+  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

%BIC selected model + cumulative success/cumulative reward main effect and
%int with eachother only
model_spec_moderate_remove_proability = 'LP_Durations_All ~total_reward_indicator_prob_All:total_reward_indicator_All +  total_reward_indicator_prob_All + total_reward_indicator_All+  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

lme__moderate_remove_probability = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability)
lme__moderate_remove_probability_coef = lme__moderate_remove_probability.Coefficients.Estimate;
lme__moderate_remove_probability_pvals = lme__moderate_remove_probability.Coefficients.pValue;	
lme__moderate_remove_probability_se = lme__moderate_remove_probability.Coefficients.SE;
lme__moderate_remove_probability_names = lme__moderate_remove_probability.Coefficients.Name;
prob_anova = anova(lme__moderate_remove_probability);
prob_anova_F = prob_anova.FStat;
prob_anova_pval = prob_anova.pValue;
[b bnames bstats] = randomEffects(lme__moderate_remove_probability);

% lme__moderate_remove_probability = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability,'DummyVarCoding','effects')

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

%BIC selected model
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%BIC selected model + cumulative success and rewar
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All*total_reward_indicator_All+ prob_indicator:total_reward_indicator_All + prob_indicator:total_reward_indicator_prob_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%remove success + reward int
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All + total_reward_indicator_All+ prob_indicator:total_reward_indicator_All + prob_indicator:total_reward_indicator_prob_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%add 3-way cumulative success/reward:n1: prob group
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~total_reward_indicator_prob_All:n_minus_one_Durations_All:prob_indicator + total_reward_indicator_All:n_minus_one_Durations_All:prob_indicator + total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All + total_reward_indicator_All+ prob_indicator:total_reward_indicator_All + prob_indicator:total_reward_indicator_prob_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%main effects and int with prob gorup only for cumulative vars
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~total_reward_indicator_prob_All:prob_indicator + total_reward_indicator_All:prob_indicator + total_reward_indicator_prob_All:total_reward_indicator_All:prob_indicator + total_reward_indicator_prob_All:total_reward_indicator_All +  total_reward_indicator_prob_All + total_reward_indicator_All+  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';

lme__moderate_remove_probability_probgroup = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_proability_probgroup,'DummyVarCoding','reference')
lme__moderate_remove_probability_probgroup_coef = lme__moderate_remove_probability_probgroup.Coefficients.Estimate;
lme__moderate_remove_probability_probgroup_pvals = lme__moderate_remove_probability_probgroup.Coefficients.pValue;	
lme__moderate_remove_probability_probgroup_se = lme__moderate_remove_probability_probgroup.Coefficients.SE;
lme__moderate_remove_probability_probgroup_names = lme__moderate_remove_probability_probgroup.Coefficients.Name;
[b bnames bstats] = randomEffects(lme__moderate_remove_probability_probgroup);

probgroups_anova = anova(lme__moderate_remove_probability_probgroup);
probgroups_anova_F = probgroups_anova.FStat;
probgroups_anova_pval = probgroups_anova.pValue;


%BIC based selection of cumulative reward stuff to see if it improves model
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
%BIC selected model + cumulative success and rewar
model_spec_moderate_remove_proability_probgroup = 'LP_Durations_All ~total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All*total_reward_indicator_All+ prob_indicator:total_reward_indicator_All + prob_indicator:total_reward_indicator_prob_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';



%BIC based selection of cumulative reward stuff to see if it improves model
model_spec_prob_baseline = 'LP_Durations_All ~prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_baseline = fitlme(T10_Logical_and_Continuous,model_spec_prob_baseline,'DummyVarCoding','reference');
%now add cumulative success me
model_spec_prob_cum_success = 'LP_Durations_All ~total_reward_indicator_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_success = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_success,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_success)
%no impovement
%now add cum rew me instead
model_spec_prob_cum_reward = 'LP_Durations_All ~total_reward_indicator_prob_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_reward = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_reward,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_reward)
%no improvement

%add both
model_spec_prob_cum_reward_success = 'LP_Durations_All ~total_reward_indicator_prob_All + total_reward_indicator_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_reward_success = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_reward_success,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_reward_success)
%no improvement

%add both + n1:success
model_spec_prob_cum_reward_success_n1 = 'LP_Durations_All ~total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All + total_reward_indicator_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_reward_success_n1 = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_reward_success_n1,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_reward_success_n1)
%no improvement in bic, though p-val sig


%add both + n1:reward
model_spec_prob_cum_reward_n1_success_n1 = 'LP_Durations_All ~total_reward_indicator_prob_All:n_minus_one_Durations_All +total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All + total_reward_indicator_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_reward_n1_success_n1 = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_reward_n1_success_n1,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_reward_n1_success_n1)
%no improvement in bic, though p-val sig

%add 3 way ints with prob
model_spec_prob_cum_reward_n1_success_n1_probint = 'LP_Durations_All ~total_reward_indicator_prob_All:n_minus_one_Durations_All:prob_indicator + total_reward_indicator_All:n_minus_one_Durations_All:prob_indicator + total_reward_indicator_prob_All:n_minus_one_Durations_All +total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All + total_reward_indicator_All + prob_indicator + prob_indicator:reward_indicator_n1_All + prob_indicator:reward_indicator_n1_All:n_minus_one_Durations_All + moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator + (1|mouse_indicator)+(1|day_indicator)';
prob_cum_reward_n1_success_n1_probint = fitlme(T10_Logical_and_Continuous,model_spec_prob_cum_reward_n1_success_n1_probint,'DummyVarCoding','reference');
compare(prob_baseline,prob_cum_reward_n1_success_n1_probint)
%no improvement in bic, though p-val sig

% Create an LM version for graphical purposes
 model_spec_moderate_remove_proability_lm = 'LP_Durations_All ~ total_reward_indicator_All:moving_average_lp_length_n7andback_All + total_reward_indicator_All:n_minus_one_Durations_All + total_reward_indicator_prob_All:moving_average_lp_length_n7andback_All + total_reward_indicator_prob_All:n_minus_one_Durations_All +  total_reward_indicator_prob_All*total_reward_indicator_All+  moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:criteria_percent_indicator +   moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:reward_indicator_n1_All   + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:criteria_percent_indicator +    n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:reward_indicator_n1_All  +  HE_n_1_Indicator_All + ipi1 + reward_indicator_n1_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All + criteria_percent_indicator';
simple_lm__moderate_remove_probability = fitlm(T10_Logical_and_Continuous,model_spec_moderate_remove_proability_lm)

figure()
plotPartialDependence(lme__moderate_remove_probability,{'n_minus_one_Durations_All','reward_indicator_n1_All'})

figure()
plotPartialDependence(lme__moderate_remove_probability,{'moving_average_lp_length_n7andback_All','reward_indicator_n1_All'})

T10_Logical_and_Continuous.reward_indicator_n1_All =double(T10_Logical_and_Continuous.reward_indicator_n1_All)
onesss = find(T10_Logical_and_Continuous.reward_indicator_n1_All ==1)

figure()
plotInteraction(simple_lm__moderate_remove_probability,'n_minus_one_Durations_All','reward_indicator_n1_All')

figure()
plotInteraction(simple_lm__moderate_remove_probability,'moving_average_lp_length_n7andback_All','reward_indicator_n1_All')
 

%% Save all the data for later use if wanted. Need to save it in old matlab format
%file format (v7.3) due to size.
% save('1000shufdata-33-3-no5-7' , '-v7.3')

toc   