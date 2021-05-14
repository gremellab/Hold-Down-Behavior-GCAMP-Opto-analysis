function [GCAMP] = GCAMP_Linear_Regressions_r(GCAMP, base_time_start, base_time_end, time_end)

%create regressions to predict GCAMP activity given:
%current press duration
%n back press durations
%n back outcomes, HE, IPI, Time, %correct

%start with relatively simple model
Durations = GCAMP.HoldDown_times;
%calculate moving avg. from n-2 to n-116
 ma_length = 60;
moving_average_lp_length = movmean(Durations,[ma_length,0]);
moving_average_lp_length = [NaN; NaN; moving_average_lp_length];
moving_average_lp_length = moving_average_lp_length(1:end-2);
%n7 and further back moving avg to exlude presses n1-n6
moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; moving_average_lp_length];
moving_average_lp_length_n7andback = moving_average_lp_length_n7andback(1:end-5);
          
        %% n-x lever press lengths
       n_back_Lengths = {};
           n_lengths_array = [];
           Logical_n_lengths_array= [];
           Logical_n_back_Lengths ={};
          press_indices = 1:length(Durations);
           press_indices = press_indices';
                      %for n-backs 1 to 10
       for x = 1:10 
       for n_it = 1:length(Durations)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           elseif  press_indices(n_it) > x
           n_back = Durations(n_it- x);               
           end
           n_lengths_array = [n_lengths_array n_back];
           
           %add in logical lengths 0/1 for failure/success
           Logical_n_lengths_array = n_lengths_array >= GCAMP.Criteria;
                       
       end
             
         n_back_Lengths = [n_back_Lengths; n_lengths_array]; 
         n_lengths_array =[];
         Logical_n_back_Lengths = [Logical_n_back_Lengths; Logical_n_lengths_array];  
         Logical_n_lengths_array =[];
       end
       
       %the logical tests turn NaNs into 0, but we need the NaNs to pad the
       %n-back array
       %loop through the array and convert the first x logical values to
       %Nans to make sure they line up with the absolute values
       for i = 1:length(Logical_n_back_Lengths)
           Logical_n_back_Lengths{i,1} =double(Logical_n_back_Lengths{i,1});
           Logical_n_back_Lengths{i,1}(1,1:i) = NaN;
       end
    
%% create data table with the behavioral info.
%can do within a single animal/day linear regression just to check stuff
%out
n_back_All = [];
for i = 1:length(n_back_Lengths)
    n_back_All = [n_back_All; n_back_Lengths{i}];
end
n_back_All = n_back_All';

  data_table = array2table(n_back_All(:,1:10),'VariableNames',{'n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All', 'n_minus_four_Durations_All',...
    'n_minus_five_Durations_All', 'n_minus_six_Durations_All', 'n_minus_seven_Durations_All',...
    'n_minus_eight_Durations_All', 'n_minus_nine_Durations_All', 'n_minus_ten_Durations_All'});

data_table.n_minus_zero_Durations_All = Durations;
%re-arrange order so n-0 is first
data_table = [data_table(:,end) data_table(:,1) data_table(:,2:end-1)];
%add in all other beh variables
%moving average, hen-1, criteria%, n-1 upstate, n-1 reward or not, hen-2,
%ipi1, ipi2
%this lines up ipi with n-1
ipi1 = [NaN; diff(GCAMP.LP_ON_timestamps)];
ipi2 = [NaN; ipi1(1:end-1)];
data_table.ipi1 = ipi1;
data_table.ipi2 = ipi2;

data_table.n_minus_one_reward_All =  data_table.n_minus_one_Durations_All >=GCAMP.Criteria;
data_table.n_minus_one_reward_All =categorical(data_table.n_minus_one_reward_All);
data_table.moving_average_lp_length_n7andback = moving_average_lp_length_n7andback;

logical_rewards = Durations >=GCAMP.Criteria;
criteria_percent = length(Durations)/sum(logical_rewards);
criteria_percent_indicator =  ones(length(Durations),1)*criteria_percent;
data_table.criteria_percent_indicator = criteria_percent_indicator;
%the he idx is tricky. Need to paste all the LP, HE and RE timestamps into a
%single array (keeping a variable to identify them) and then sort based on
%time
LP_start_times = GCAMP.LP_ON_timestamps;
LP_start_times =[ LP_start_times repmat(10,[length(LP_start_times) 1])  ];
 
HE_start_times = GCAMP.HE_ON_timestamps;
HE_start_times =[ HE_start_times repmat(12,[length(HE_start_times) 1])  ];

RE_start_times = GCAMP.RE_ON_timestamps;
RE_start_times =[ RE_start_times repmat(20,[length(RE_start_times) 1])  ];

Event_start_times = [LP_start_times; HE_start_times; RE_start_times];
[spot og_idx ] = sort(Event_start_times(:,1));
sorted_Events = Event_start_times(og_idx,:);
GCAMP.sorted_Events = sorted_Events;
%now loop through to find when a lp 
nanpad = [NaN NaN; NaN NaN];
sorted_Events = [sorted_Events; nanpad];
HE_Indicator =[];
HE_counter_rewarded = 0;
HE_counter_general =0;
HE_counter_unrewarded =0;
for he_it = 1:length(sorted_Events)      
        
if sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ==20  %if lp and rewarded
    
    if sorted_Events(he_it+2,2) == 12 %if theres a 12 after the 20 add it to the general he counter and reward specific
        HE_counter_rewarded = HE_counter_rewarded +1;
        HE_counter_general = HE_counter_general +1;
        
        HE_Indicator = [HE_Indicator; 1];
        
    elseif sorted_Events(he_it+2,2) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
end

if sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ~=20  %if LP was unrewarded
    if sorted_Events(he_it+1,2) == 12 %if theres a 12 after the 15 add it to the general he counter and unreward specific
        HE_counter_unrewarded = HE_counter_unrewarded +1;
        HE_counter_general = HE_counter_general +1;
        
        HE_Indicator = [HE_Indicator; 1];
        
    elseif sorted_Events(he_it+1,2) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
             
end
end
HE_n1_Indicator = [NaN; HE_Indicator];
HE_n1_Indicator(end) =[];

data_table.HE_indicator = HE_Indicator;
data_table.HE_n1_Indicator = HE_n1_Indicator;
data_table.HE_indicator = categorical(data_table.HE_indicator);
data_table.HE_n1_Indicator = categorical(data_table.HE_n1_Indicator);
%we have lp timestamps, but they are in bonsai time. Basically just need to
%make them in relation to the session, so timestamp 1 will become 0, and
%then timed up from there
lp_on_times  = GCAMP.LP_ON_timestamps;
lp_on_times = lp_on_times- lp_on_times(1);
data_table.lp_on_times = lp_on_times;

%% Code to shuffle the order x number of times in order to build a distribution of shuffled data.
%also need to build the n-back arrays from these shuffled datasets
Shuffled_Durs_Distribtuion = [];
shuffled_ipi1_Distribution = [];
shuffled_ipi2_Distribution =[];
shuffled_he1_Distribution =[];
Shuffled_ts_Distribution =[];
shuffled_MA_n7_Distribution =[];
shuf_n_back_Lengths_distribution = {};
shuf_Logical_n_back_Lengths_distribution = {};

for shuf_it = 1:1000
    Shuffled_Durs = Durations;
    Shuffled_Durs = Shuffled_Durs(randperm(length(Shuffled_Durs)));
    Shuffled_Durs_Distribtuion = [Shuffled_Durs_Distribtuion Shuffled_Durs]; 
          
    Shuffled_ts = lp_on_times;
    Shuffled_ts = Shuffled_ts(randperm(length(Shuffled_ts)));
    Shuffled_ts_Distribution = [Shuffled_ts_Distribution Shuffled_ts];
       
    %calculate a moving average for each shuffle
    shuffled_moving_average_lp_length = movmean(Shuffled_Durs,[ma_length,0]);      
    shuffled_moving_average_lp_length = [NaN; NaN; shuffled_moving_average_lp_length];  
    shuffled_moving_average_lp_length = shuffled_moving_average_lp_length(1:end-2);
     
        shuffled_moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; shuffled_moving_average_lp_length];
    shuffled_moving_average_lp_length_n7andback = shuffled_moving_average_lp_length_n7andback(1:end-5);
     shuffled_MA_n7_Distribution = [shuffled_MA_n7_Distribution shuffled_moving_average_lp_length_n7andback];
  
    %independentyly shuffle ipi and he
    shuf_ipi1= ipi1(2:end);
    shuffled_ipi1 = shuf_ipi1(randperm(length(shuf_ipi1)));
    %removed the original nan, now add it back in
    shuffled_ipi1 = [NaN; shuffled_ipi1];
    shuffled_ipi1 = shuffled_ipi1(1:end-1);
    shuffled_ipi1_Distribution=[shuffled_ipi1_Distribution  shuffled_ipi1];  
     
    shuf_ipi2 = ipi2(3:end);
    shuffled_ipi2 =shuf_ipi2(randperm(length(shuf_ipi2)));
    shuffled_ipi2 = [NaN; NaN; shuffled_ipi2];
    shuffled_ipi2_Distribution=[shuffled_ipi2_Distribution  shuffled_ipi2];  

    shuf_he1 = HE_n1_Indicator(2:end);
    shuf_he1 =  shuf_he1(randperm(length(shuf_he1)));
    shuf_he1 = [NaN; shuf_he1];
    shuffled_he1_Distribution =[shuffled_he1_Distribution shuf_he1];

end

GCAMP.shuffled_ipi1_Distribution = shuffled_ipi1_Distribution;
GCAMP.shuffled_ipi2_Distribution = shuffled_ipi2_Distribution;
GCAMP.shuffled_he1_Distribution =shuffled_he1_Distribution;
GCAMP.Shuffled_ts_Distribution =Shuffled_ts_Distribution;
GCAMP.shuffled_MA_n7_Distribution =shuffled_MA_n7_Distribution;
 
%now shuffle the nback lengths  
for shuffle_length = 1:size(Shuffled_Durs_Distribtuion,2)
    shuf_n_back_Lengths = [];
    shuf_n_lengths_array = [];
    shuf_Logical_n_lengths_array= [];
    shuf_Logical_n_back_Lengths =[];
   
           current_shuffle_lps = Shuffled_Durs_Distribtuion(:,shuffle_length);
           shuf_press_indices = 1:length(current_shuffle_lps);
           shuf_press_indices = shuf_press_indices';
           
           %for n-backs 1 to 10
       for x = 1:10 
   
       for n_it = 1:length(current_shuffle_lps)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           elseif  press_indices(n_it) > x
           n_back = current_shuffle_lps(n_it- x);               
           end
           shuf_n_lengths_array = [shuf_n_lengths_array n_back];
           
           %add in logical lengths 0/1 for failure/success
          
           shuf_Logical_n_lengths_array = shuf_n_lengths_array >= GCAMP.Criteria;
           
           end
       % shuf_n_lengths_array
      shuf_n_back_Lengths = [shuf_n_back_Lengths; shuf_n_lengths_array]; 
      shuf_n_lengths_array =[];
      
      %the logical test converts all NaNs to 0, need to put them back so we
      %don't inappropriately put lever press prior to the start of the
      %session
               shuf_Logical_n_lengths_array = double(shuf_Logical_n_lengths_array);
        shuf_Logical_n_lengths_array(1:x) = NaN;
         shuf_Logical_n_back_Lengths = [shuf_Logical_n_back_Lengths; shuf_Logical_n_lengths_array];  
         shuf_Logical_n_lengths_array =[];
   
       end
       shuf_n_back_Lengths_distribution = [shuf_n_back_Lengths_distribution; shuf_n_back_Lengths];
       shuf_Logical_n_back_Lengths_distribution = [shuf_Logical_n_back_Lengths_distribution; shuf_Logical_n_back_Lengths];    
       
end
% GCAMP.shuf_n1_reward_Distribution = shuf_Logical_n_back_Lengths_distribution
GCAMP.shuf_n_back_Lengths_distribution = shuf_n_back_Lengths_distribution;
GCAMP.Shuffled_Durs_Distribtuion = Shuffled_Durs_Distribtuion;
GCAMP.shuf_Logical_n_back_Lengths_distribution =shuf_Logical_n_back_Lengths_distribution;

%% get gcamp data for relevant events. Segment into timepoints of interest and mean and AUC
%see if baseline, baseline_raw (no running filter), or z-scored raw data
%affects the direction/magnitude of any effects.

%baseline z-scored data
Onset_baselinezscore_Raw = GCAMP.baseline_z_score_LP_On_Raw;
Offset_baselinezscore_Raw = GCAMP.baseline_z_score_LP_OFF_Raw;

%lp offset for rewarded lp is the same as reward
RE_baselinezscore_raw = GCAMP.baseline_z_score_LP_OFF_Raw_Met;

%get the baselines themselves
baselines_on_raw = GCAMP.baselines_on_raw;

%now z-score the raw signal within mouse, zscore acts on columns so flip it to z-score within a lp
Onset_raw = GCAMP.raw_F_LP_ON;
Offset_raw = GCAMP.raw_F_LP_OFF;
RE_raw = GCAMP.raw_F_RE_ON;
Onset_raw = Onset_raw';
Onset_raw_zscored = GCAMP.raw_F_LP_ON_Z;

Offset_raw = Offset_raw';
Offset_raw_zscored = zscore(Offset_raw);
Offset_raw_zscored = Offset_raw_zscored';

RE_raw = RE_raw';
RE_raw_zscored = zscore(RE_raw);
RE_raw_zscored = RE_raw_zscored';

%interpolated baseline-zscored of raw signal data
interpolated_base = GCAMP.interp_data;

%now we want to take several periods around press onset and get both AVG
%and AUC
event_onset_idx = 1+abs(GCAMP.SR*base_time_end);

%% baseline raw z-scored 
%try a few different windows (.5s, 1s) and before/after onset
base_neg_1_to_onset_raw = Onset_baselinezscore_Raw(:,event_onset_idx-GCAMP.SR*1:event_onset_idx);
base_neg_point5_to_onset_raw = Onset_baselinezscore_Raw(:,event_onset_idx-GCAMP.SR*.5:event_onset_idx);
base_onset_to_pont5_raw = Onset_baselinezscore_Raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*.5);
base_onset_to1_raw = Onset_baselinezscore_Raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*1);

%% baseline onset raw
%try mean, auc, slope
base_neg_1_to_onset_AUC_raw = trapz(base_neg_1_to_onset_raw,2);
base_neg_1_to_onset_Mean_raw = mean(base_neg_1_to_onset_raw,2);
base_neg_point5_to_onset_AUC_raw = trapz(base_neg_point5_to_onset_raw,2);
base_neg_point5_to_onset_Mean_raw = mean(base_neg_point5_to_onset_raw,2);
base_onset_to_pont5_AUC_raw = trapz(base_onset_to_pont5_raw,2);
base_onset_to_pont5_mean_raw = mean(base_onset_to_pont5_raw,2);
base_onset_to1_AUC_raw = trapz(base_onset_to1_raw,2);
base_onset_to1_Mean_raw = mean(base_onset_to1_raw,2);
base_neg_1_to_onset_Slope_raw= (mean(base_neg_1_to_onset_raw(:,1:4),2)-mean(base_neg_1_to_onset_raw(:,end-4:end),2))/ (17);
base_onset_to1_Slope_raw= (mean(base_onset_to1_raw(:,1:4),2)-mean(base_onset_to1_raw(:,end-4:end),2))/ (17);
base_neg_1_to_onset_std_raw = std(base_neg_1_to_onset_raw,0,2);
data_table.base_neg_1_to_onset_std_raw = base_neg_1_to_onset_std_raw ;

data_table.base_neg_1_to_onset_AUC_raw = base_neg_1_to_onset_AUC_raw;
data_table.base_neg_1_to_onset_Mean_raw = base_neg_1_to_onset_Mean_raw;
data_table.base_neg_point5_to_onset_AUC_raw = base_neg_point5_to_onset_AUC_raw;
data_table.base_neg_point5_to_onset_Mean_raw =base_neg_point5_to_onset_Mean_raw;
data_table.base_onset_to_pont5_AUC_raw = base_onset_to_pont5_AUC_raw;
data_table.base_onset_to_pont5_mean_raw = base_onset_to_pont5_mean_raw;
data_table.base_onset_to1_AUC_raw = base_onset_to1_AUC_raw;
data_table.base_onset_to1_Mean_raw = base_onset_to1_Mean_raw;
data_table.braw_neg_1_to_onset_Slope_raw= base_neg_1_to_onset_Slope_raw;
data_table.raw_onset_to1_Slope_raw= base_onset_to1_Slope_raw;
%% lp offset
%now offset. We see sig differences at about 2-5s post offset (reward), so
%add that in as well as 0 to 1s
%baseline raw
base_neg_1_to_offset_raw = Offset_baselinezscore_Raw(:,event_onset_idx-GCAMP.SR*1:event_onset_idx);
base_neg_point5_to_offset_raw = Offset_baselinezscore_Raw(:,event_onset_idx-GCAMP.SR*.5:event_onset_idx);
base_offset_to_pont5_raw = Offset_baselinezscore_Raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*.5);
base_offset_to1_raw = Offset_baselinezscore_Raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*1);
base_offset_2_to_5_raw = Offset_baselinezscore_Raw(:,event_onset_idx+GCAMP.SR*2:event_onset_idx+GCAMP.SR*5);


%baseline offset raw
base_neg_1_to_offset_AUC_raw = trapz(base_neg_1_to_offset_raw,2);
base_neg_1_to_offset_Mean_raw = mean(base_neg_1_to_offset_raw,2);
base_neg_point5_to_offset_AUC_raw = trapz(base_neg_point5_to_offset_raw,2);
base_neg_point5_to_offset_Mean_raw = mean(base_neg_point5_to_offset_raw,2);
base_offset_to_pont5_AUC_raw = trapz(base_offset_to_pont5_raw,2);
base_offset_to_pont5_mean_raw = mean(base_offset_to_pont5_raw,2);
base_offset_to1_AUC_raw = trapz(base_offset_to1_raw,2);
base_offset_to1_Mean_raw = mean(base_offset_to1_raw,2);
base_offset_to1_std_raw = std(base_offset_to1_raw,0,2);

base_offset_2_to_5_AUC_raw = trapz(base_offset_2_to_5_raw,2);
base_offset_2_to_5_Mean_raw = mean(base_offset_2_to_5_raw,2);
base_neg_1_to_offset_Slope_raw= (mean(base_neg_1_to_offset_raw(:,1:4),2)-mean(base_neg_1_to_offset_raw(:,end-4:end),2))/ (17);
base_offset_to1_Slope_raw= (mean(base_offset_to1_raw(:,1:4),2)-mean(base_offset_to1_raw(:,end-4:end),2))/ (17);
base_offset_2_to_5_Slope_raw= (mean(base_offset_2_to_5_raw(:,1:4),2)-mean(base_offset_2_to_5_raw(:,end-4:end),2))/ (17);

data_table.base_offset_to1_std_raw = base_offset_to1_std_raw;
data_table.base_neg_1_to_offset_AUC_raw =base_neg_1_to_offset_AUC_raw;
data_table.base_neg_1_to_offset_Mean_raw =base_neg_1_to_offset_Mean_raw;
data_table.base_neg_point5_to_offset_AUC_raw = base_neg_point5_to_offset_AUC_raw;
data_table.base_neg_point5_to_offset_Mean_raw = base_neg_point5_to_offset_Mean_raw;
data_table.base_offset_to_pont5_AUC_raw = base_offset_to_pont5_AUC_raw;
data_table.base_offset_to_pont5_mean_raw = base_offset_to_pont5_mean_raw;
data_table.base_offset_to1_AUC_raw = base_offset_to1_AUC_raw;
data_table.base_offset_to1_Mean_raw = base_offset_to1_Mean_raw;
data_table.base_offset_2_to_5_AUC_raw = base_offset_2_to_5_AUC_raw;
data_table.base_offset_2_to_5_Mean_raw = base_offset_2_to_5_Mean_raw;
data_table.base_offset_2_to_5_Slope_raw = base_offset_2_to_5_Slope_raw;
data_table.base_offset_to1_Slope_raw = base_offset_to1_Slope_raw;
data_table.base_neg_1_to_offset_Slope_raw = base_neg_1_to_offset_Slope_raw;

%% reward on/offset
%baseline raw
base_RE_raw_neg_1_to_onset = RE_baselinezscore_raw(:,event_onset_idx-GCAMP.SR*1:event_onset_idx);
base_RE_raw_neg_point5_to_onset = RE_baselinezscore_raw(:,event_onset_idx-GCAMP.SR*.5:event_onset_idx);
base_RE_raw_onset_to_pont5 = RE_baselinezscore_raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*.5);
base_RE_raw_onset_to1 = RE_baselinezscore_raw(:,event_onset_idx:event_onset_idx+GCAMP.SR*1);

%mean and AUC base raw
base_RE_neg_1_to_offset_AUC_raw = trapz(base_RE_raw_neg_1_to_onset,2);
base_RE_neg_1_to_offset_Mean_raw = mean(base_RE_raw_neg_1_to_onset,2);
base_RE_neg_point5_to_offset_AUC_raw = trapz(base_RE_raw_neg_point5_to_onset,2);
base_RE_neg_point5_to_offset_Mean_raw = mean(base_RE_raw_neg_point5_to_onset,2);
base_RE_offset_to_pont5_AUC_raw = trapz(base_RE_raw_onset_to_pont5,2);
base_RE_offset_to_pont5_mean_raw = mean(base_RE_raw_onset_to_pont5,2);
base_RE_offset_to1_AUC_raw = trapz(base_RE_raw_onset_to1,2);
base_RE_offset_to1_Mean_raw = mean(base_RE_raw_onset_to1,2);

%% interpolated presses
%divide the interpolated presses into 4 bins 
interp_quarter = size(interpolated_base,2)/4;
interp_first_quart = interpolated_base(:,1:interp_quarter);
interp_second_quart = interpolated_base(:,interp_quarter+1:interp_quarter*2);
interp_third_quart = interpolated_base(:,1+interp_quarter*2:interp_quarter*3);
interp_fourth_quart = interpolated_base(:,1+interp_quarter*3:end);

%get mean and AUC
interp_all_AUC = trapz(interpolated_base,2);
interp_all_mean =mean(interpolated_base,2);
interp_firstquart_AUC = trapz(interp_first_quart,2);
interp_firstquart_mean =mean(interp_first_quart,2);
interp_secondquart_AUC = trapz(interp_second_quart,2);
interp_secondquart_mean =mean(interp_second_quart,2);
interp_thirdquart_AUC = trapz(interp_third_quart,2);
interp_thirdquart_mean =mean(interp_third_quart,2);
interp_fourthquart_AUC = trapz(interp_fourth_quart,2);
interp_fourthquart_mean =mean(interp_fourth_quart,2);
interp_all_Slope = (nanmean(interpolated_base(:,1:4),2)-nanmean(interpolated_base(:,end-4:end),2))/ (16);
interp_all_std = std(interpolated_base,0,2);

data_table.interp_all_std = interp_all_std;
data_table.interp_all_AUC = interp_all_AUC;
data_table.interp_all_mean =interp_all_mean;
data_table.interp_firstquart_AUC = interp_firstquart_AUC;
data_table.interp_firstquart_mean =interp_firstquart_mean;
data_table.interp_secondquart_AUC = interp_secondquart_AUC;
data_table.interp_secondquart_mean =interp_secondquart_mean;
data_table.interp_thirdquart_AUC = interp_thirdquart_AUC;
data_table.interp_thirdquart_mean =interp_thirdquart_mean;
data_table.interp_fourthquart_AUC = interp_fourthquart_AUC;
data_table.interp_fourthquart_mean =interp_fourthquart_mean;
data_table.interp_all_Slope = interp_all_Slope;

%% baselines

baselines_raw_mean = mean(baselines_on_raw,2);
baselines_raw_AUC = trapz(baselines_on_raw,2);
baselines_raw_slope = (nanmean(baselines_on_raw(:,1:4),2)-nanmean(baselines_on_raw(:,end-4:end),2))/ (197);

data_table.baselines_raw_mean = baselines_raw_mean;
data_table.baselines_raw_AUC = baselines_raw_AUC;
data_table.baselines_raw_slope = baselines_raw_slope;

%% add in lagged variables of the chief variables of interest
data_table.interp_all_AUC_n1 = [NaN; data_table.interp_all_AUC(1:end-1)]; 
data_table.interp_all_AUC_n2 = [NaN; data_table.interp_all_AUC_n1(1:end-1)]; 
data_table.interp_all_AUC_n3 = [NaN; data_table.interp_all_AUC_n2(1:end-1)]; 
data_table.interp_all_AUC_n4 = [NaN; data_table.interp_all_AUC_n3(1:end-1)]; 
data_table.interp_all_AUC_n5 = [NaN; data_table.interp_all_AUC_n4(1:end-1)]; 
data_table.interp_all_AUC_n6 = [NaN; data_table.interp_all_AUC_n5(1:end-1)]; 

data_table.base_neg_1_to_onset_AUC_raw_n1 = [NaN; data_table.base_neg_1_to_onset_AUC_raw(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n2 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n1(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n3 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n2(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n4 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n3(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n5 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n4(1:end-1)]; 
data_table.base_neg_1_to_onset_AUC_raw_n6 = [NaN; data_table.base_neg_1_to_onset_AUC_raw_n5(1:end-1)]; 

data_table.base_offset_to1_AUC_raw_n1 = [NaN; data_table.base_offset_to1_AUC_raw(1:end-1)];
data_table.base_offset_to1_AUC_raw_n2 = [NaN; data_table.base_offset_to1_AUC_raw_n1(1:end-1)];
data_table.base_offset_to1_AUC_raw_n3 = [NaN; data_table.base_offset_to1_AUC_raw_n2(1:end-1)];
data_table.base_offset_to1_AUC_raw_n4 = [NaN; data_table.base_offset_to1_AUC_raw_n3(1:end-1)];
data_table.base_offset_to1_AUC_raw_n5 = [NaN; data_table.base_offset_to1_AUC_raw_n4(1:end-1)];
data_table.base_offset_to1_AUC_raw_n6 = [NaN; data_table.base_offset_to1_AUC_raw_n5(1:end-1)];

data_table.base_offset_2_to_5_AUC_raw_n1 = [NaN; data_table.base_offset_2_to_5_AUC_raw(1:end-1)];
data_table.base_offset_2_to_5_AUC_raw_n2 = [NaN; data_table.base_offset_2_to_5_AUC_raw_n1(1:end-1)];
data_table.base_offset_2_to_5_AUC_raw_n3 = [NaN; data_table.base_offset_2_to_5_AUC_raw_n2(1:end-1)];
data_table.base_offset_2_to_5_AUC_raw_n4 = [NaN; data_table.base_offset_2_to_5_AUC_raw_n3(1:end-1)];
data_table.base_offset_2_to_5_AUC_raw_n5 = [NaN; data_table.base_offset_2_to_5_AUC_raw_n4(1:end-1)];
data_table.base_offset_2_to_5_AUC_raw_n6 = [NaN; data_table.base_offset_2_to_5_AUC_raw_n5(1:end-1)];
%%
%now we need to save stuff into the GCAMP structure
%in particular need to save the table structures, and add in indicator
%variables for mouse and day of training
reps = size(data_table,1);

mouse_indicator = repmat(GCAMP.onlyID ,reps,1);
day_indicator = repmat(GCAMP.training_day,reps,1);
criteria_indicator = repmat(GCAMP.Criteria,reps,1);
data_table.mouse_indicator = mouse_indicator;
data_table.day_indicator = day_indicator;
data_table.criteria_indicator = criteria_indicator;


data_table_variables = data_table.Properties.VariableNames;
GCAMP.data_table_variables = data_table_variables;
%need to convert back to an array so that we can past all thetables
%together 
regression_cell= table2cell(data_table);
GCAMP.regression_cell = regression_cell;

%% 
%Now fit individual mouse data
longidx = data_table.n_minus_zero_Durations_All>=10000;
%need to add n0 var to table
data_table.n_zero_reward = data_table.n_minus_zero_Durations_All>=data_table.criteria_indicator;

%simple beh with just durations
modelspec_behavior_lponly = 'n_minus_zero_Durations_All ~ n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All + n_minus_six_Durations_All + n_minus_seven_Durations_All + n_minus_eight_Durations_All + n_minus_nine_Durations_All + n_minus_ten_Durations_All';
lme_behavior_lponly = fitlme(data_table,modelspec_behavior_lponly,'Exclude',longidx);
%33-3 behavioral model full
modelspec_behavior_reduced_333 = 'n_minus_zero_Durations_All ~n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2+n_minus_two_Durations_All*ipi2+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+ lp_on_times:moving_average_lp_length_n7andback';
lme_behavior_reduced_333 = fitlme(data_table,modelspec_behavior_reduced_333,'Exclude',longidx);

%activity duration simple interpolated duration
modelspec_interp_all_n6_act = 'interp_all_AUC ~ n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All';
lme_interp_all_n6_act = fitlme(data_table,modelspec_interp_all_n6_act,'Exclude',longidx);
%activity full beh model
modelspec_interp_all_n6_act_andints = 'interp_all_AUC ~ n_minus_zero_Durations_All + n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+ lp_on_times:moving_average_lp_length_n7andback';
lme_interp_all_n6_act_andints = fitlme(data_table,modelspec_interp_all_n6_act_andints,'Exclude',longidx);

%Pre onset simple
modelspec_pre_onset_base_raw_1s_n6_act = 'base_neg_1_to_onset_AUC_raw ~  n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All ';
lme_pre_onset_base_raw_1s_n6_act = fitlme(data_table,modelspec_pre_onset_base_raw_1s_n6_act,'Exclude',longidx);
%pre onset complex
modelspec_preonset_all_n6_act_andints = 'base_neg_1_to_onset_AUC_raw ~ n_minus_zero_Durations_All + n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All + lp_on_times:moving_average_lp_length_n7andback';
lme_preonset_all_n6_act_andints = fitlme(data_table,modelspec_preonset_all_n6_act_andints,'Exclude',longidx);

%post offset simple
modelspec_post_offset_base_raw_to1_n6_act = 'base_offset_to1_AUC_raw ~  n_minus_zero_Durations_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All';
lme_post_offset_base_raw_to1_n6_act = fitlme(data_table,modelspec_post_offset_base_raw_to1_n6_act,'Exclude',longidx);
%post offset complex
modelspec_postoffset_all_n6_act_andints_n0 = 'base_offset_to1_AUC_raw ~n_zero_reward:moving_average_lp_length_n7andback + n_minus_zero_Durations_All * n_zero_reward  + n_minus_zero_Durations_All + n_minus_six_Durations_All+n_minus_five_Durations_All+n_minus_four_Durations_All+n_minus_three_Durations_All+ipi2*n_minus_two_Durations_All+n_minus_one_reward_All*n_minus_one_Durations_All+ n_minus_one_Durations_All*HE_n1_Indicator+lp_on_times*n_minus_one_Durations_All+ipi1*n_minus_one_Durations_All+moving_average_lp_length_n7andback + moving_average_lp_length_n7andback:HE_n1_Indicator + ipi1:moving_average_lp_length_n7andback +  moving_average_lp_length_n7andback:n_minus_one_reward_All+ lp_on_times:moving_average_lp_length_n7andback';
lme_postoffset_all_n6_act_andints_n0 = fitlme(data_table,modelspec_postoffset_all_n6_act_andints_n0,'Exclude',longidx);

GCAMP.simple_beh_r2 = lme_behavior_lponly.Rsquared.Ordinary;
GCAMP.complex_beh_r2 = lme_behavior_reduced_333.Rsquared.Ordinary;
GCAMP.simple_interp_r2 = lme_interp_all_n6_act.Rsquared.Ordinary;
GCAMP.complex_interp_r2 = lme_interp_all_n6_act_andints.Rsquared.Ordinary;
GCAMP.simple_pre_r2 = lme_pre_onset_base_raw_1s_n6_act.Rsquared.Ordinary;
GCAMP.complex_pre_r2 = lme_preonset_all_n6_act_andints.Rsquared.Ordinary;
GCAMP.simple_post_r2 = lme_post_offset_base_raw_to1_n6_act.Rsquared.Ordinary;
GCAMP.complex_post_r2 = lme_postoffset_all_n6_act_andints_n0.Rsquared.Ordinary;
GCAMP.criteria_prop = data_table.criteria_percent_indicator;
end


