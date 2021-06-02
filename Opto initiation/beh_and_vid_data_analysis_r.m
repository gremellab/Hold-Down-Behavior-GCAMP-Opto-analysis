function [Vid_Data] = beh_and_vid_data_analysis_r(beh_data_raw, vid_data_raw, opto_data_raw,zone_data_raw, Criteria, mouseID, training_day)
%last updated 5/12/21 Drew Schreienr
%Take the imported data and analyze it
%looking for things like LPs, HEs, timeing, stimulation

onlyID =mouseID(1:4);
day_indicator = mouseID(12:13);
%% session only
%there is a session on indicator in the behavioral data, so we want to find
%that and use it to truncate the other data
%code to only get when the session is on, indicated by 1 in column 1 of
% indices = find(beh_data_raw(:,1)==0);
%create new versions
beh_data = beh_data_raw;
vid_data = vid_data_raw;
opto_data = opto_data_raw;
zone_data = zone_data_raw;

%the 2nd column is the session indicator, truncate all data to only include
%from start to end of session.
start_session_idx = find(beh_data(:,2) == 1,1);
end_session_idx = find(beh_data(:,2) == 1,1,'last');

beh_data = beh_data(start_session_idx:end_session_idx,:);
start_session_ts = beh_data(1,1);
end_session_ts = beh_data(end,1);
 
%we can use the idx for the beh data, but need to use the timestamps for
%the other csvs
start_idx_opto = find(opto_data(:,1) >= start_session_ts,1);
end_idx_opto =   find(opto_data(:,1) <= end_session_ts,1,'last');
opto_data = opto_data(start_idx_opto:end_idx_opto,:);

start_idx_zone = find(zone_data(:,1) >= start_session_ts,1);
end_idx_zone = find(zone_data(:,1) <= end_session_ts,1,'last');
zone_data = zone_data(start_idx_zone:end_idx_zone,:);

start_idx_track = find(vid_data(:,1) >= start_session_ts,1);
end_idx_track = find(vid_data(:,1) <= end_session_ts,1,'last');
vid_data = vid_data(start_idx_track:end_idx_track,:);

%%
%get all the relevant behavioral variables (LP,HE,Reward)

%column 1 = timestamps
%col 2 = Session indicator (houselight)
%col 3 = LPs
%col 4 = stimulation
%col 5 = HE
% Subtract logicals from each other in Lever Press column
switches_LP = diff(beh_data(:,3));
% Subtract logicals from each other in Reinforcement column
Switches_stim = diff(beh_data(:,4));
% Subtract logicals from each other in HeadEntries column
switches_HE = diff(beh_data(:,5));
% When the difference is 1, that's event ON
LP_ON_index = find(switches_LP == 1)+1;
Stim_ON_index = find(Switches_stim == 1)+1;
HE_ON_index = find(switches_HE == 1)+1;
% When the difference is -1, that's event OFF
LP_OFF_index = find(switches_LP == -1)+1;
Stim_OFF_index = find(Switches_stim == -1)+1;
HE_OFF_index = find(switches_HE == -1)+1;
% Get the time stamps from those event indices
LP_ON_timestamps = beh_data(LP_ON_index, 1);
LP_OFF_timestamps = beh_data(LP_OFF_index, 1);
Stim_ON_timestamps = beh_data(Stim_ON_index, 1);
Stim_OFF_timestamps = beh_data(Stim_OFF_index, 1);
HE_ON_timestamps = beh_data(HE_ON_index, 1);
HE_OFF_timestamps = beh_data(HE_OFF_index,1);

%Want Matrix of times where Column 1 is event ON, Column 2 is event OFF
%sometimes there is an uneven on/off, truncate the final one if so
if length(LP_ON_timestamps) > length(LP_OFF_timestamps)
    LP_ON_timestamps(end) = [];
    LP_ON_index(end) =[];
end

LP_All = [LP_ON_timestamps LP_OFF_timestamps];
% Hold down time for each lever press
HD_time=(LP_All(:,2)-LP_All(:,1))-20;
HD_time_unrounded = HD_time;
%round hold down time to the nearest 10ms since this is med-pc's resolution
HD_time = round(HD_time,-1);

%exclude durations that are 0/neg
exclude_idx = find(HD_time < 10) ;

LP_ON_index(exclude_idx) =[];
LP_OFF_index(exclude_idx) =[];
LP_ON_timestamps(exclude_idx) =[];
LP_OFF_timestamps(exclude_idx) =[];
HD_time(exclude_idx) =[];
LP_All(exclude_idx) =[];
% Lever presses that met criteria and were reinforced
Criteria_met = find(HD_time >=  Criteria); %Criteria);
Criteria_fail = find(HD_time <   Criteria);
Total_Reinforcers = length(Criteria_met);

%% what we most care about is finding how stim affects:
%latency to make a lp
%prob of making a lp 
%duration of a lp
%similarity to preceding lever press (LME betas)

%first to find latency, we find when 0s (no stim) or 1s (yes stim) occur in
%the opto_data. Then take that timestamp and find the nearest LP timestamp
%in the analog data
opto_stim = opto_data(:,2) ==1;
opto_stim_ts = opto_data(opto_stim);
opto_no_stim = opto_data(:,2) ==0;
opto_no_stim_ts = opto_data(opto_no_stim);

[nearest_stim_idx nearest_stim_latency] =nearestpoint(opto_stim_ts(:,1),LP_ON_timestamps,'next');
[nearest_no_stim_idx nearest_no_stim_latency] =nearestpoint(opto_no_stim_ts(:,1),LP_ON_timestamps,'next');

mean_nearest_stim_latency = nanmean(nearest_stim_latency);
mean_nearest_no_stim_latency = nanmean(nearest_no_stim_latency);

%problem: there are times when animals enter zone and get or dont get stim
%multiple times before lping. We need to make sure we only take the MOST
%recent zone entrance (whether stim or not) when the same index exists.
[nearest_lp_idx nearest_lp_latency] =nearestpoint(opto_data(:,1),LP_ON_timestamps,'next');
nearest_lp_idx = [nearest_lp_idx opto_data(:,2)];
%ia is the index of the final lp before a lp
[C,ia,ic] = unique(nearest_lp_idx(:,1),'last');
nearest_lp_idx_first = nearest_lp_idx(ia,:);
nearest_lp_latency_first = nearest_lp_latency(ia);
nearest_lp_latency_first=rmmissing(nearest_lp_latency_first);

%% optional code to only look at lps that occur within x amount of time. This can be
%used to segment out data for only when there is or is not overlap with the
%LP itself 
% nearest_lp_latency_first_shortidx = nearest_lp_latency_first < 1000 ;
% nearest_lp_idx_first = nearest_lp_idx_first(nearest_lp_latency_first_shortidx,:);
% nearest_lp_latency_first = nearest_lp_latency_first(nearest_lp_latency_first_shortidx);
%%
nearest_lp_idx_first = rmmissing(nearest_lp_idx_first);

%Segment out stim/no stim latency and get means
nearest_lp_latency_first_stim = nearest_lp_latency_first(nearest_lp_idx_first(:,2)==1,:);
nearest_lp_latency_first_no_stim = nearest_lp_latency_first(nearest_lp_idx_first(:,2)==0,:);
nearest_lp_latency_first_stim_mean = nanmean(nearest_lp_latency_first_stim);
nearest_lp_latency_first_no_stim_mean = nanmean(nearest_lp_latency_first_no_stim);

%segment into stim or no bassed on 1 or 0 in col 2
nearest_lp_idx_first_stim = nearest_lp_idx_first(nearest_lp_idx_first(:,2)==1,:);
nearest_lp_idx_first_no_stim = nearest_lp_idx_first(nearest_lp_idx_first(:,2)==0,:);

%This will be an indicator variable where 1 = no stim preceding a lp (but
%still, animals met the criteria to receive stim if it weren't
%probabalistic, thus we make a direct comparison between 1st presses upon
%zone entrance without stim (1) and with stim (2).
stim_and_zone_indicator = NaN(length(HD_time),1);
stim_and_zone_indicator(nearest_lp_idx_first_no_stim(:,1)) = 1;
stim_and_zone_indicator(nearest_lp_idx_first_stim(:,1)) = 2;

%Get latencies between stimulation and lever pressing 
stim_and_zone_latency = NaN(length(HD_time),1);
nearest_lp_idx_and_latency = [nearest_lp_idx_first(:,1) nearest_lp_latency_first];
stim_and_zone_latency(nearest_lp_idx_and_latency(:,1)) =nearest_lp_idx_and_latency(:,2);

%need to match this to lp durations
%because there are oftentimes many lps made aftera a zone entrace, just get
%the unique idx to get the first lp
nearest_stim_idx_single = unique(nearest_stim_idx,'last');
nearest_no_stim_idx_single = unique(nearest_no_stim_idx,'last');
nearest_stim_idx_single = rmmissing(nearest_stim_idx_single);
nearest_no_stim_idx_single = rmmissing(nearest_no_stim_idx_single);

stim_Durations = HD_time(stim_and_zone_indicator==2);
no_stim_Durations = HD_time(stim_and_zone_indicator==1);
mean_stim_Durations=mean(stim_Durations);
mean_no_stim_Durations = mean(no_stim_Durations);

%save
Vid_Data.nearest_lp_latency_first_stim_mean = nearest_lp_latency_first_stim_mean;
Vid_Data.nearest_lp_latency_first_no_stim_mean = nearest_lp_latency_first_no_stim_mean;
Vid_Data.nearest_lp_latency_first_stim = nearest_lp_latency_first_stim;
Vid_Data.nearest_lp_latency_first_no_stim = nearest_lp_latency_first_no_stim;
Vid_Data.stim_Durations = stim_Durations;
Vid_Data.no_stim_Durations = no_stim_Durations;
Vid_Data.mean_stim_Durations = mean_stim_Durations;
Vid_Data.mean_no_stim_Durations = mean_no_stim_Durations;
Vid_Data.stim_and_zone_indicator =stim_and_zone_indicator;
Vid_Data.stim_and_zone_latency = stim_and_zone_latency;

%%
%Rather than just looking at nearest LP to stim (or no stim), we can also
%segment it based on their being a LP within the zone prior to zone leaving
current_lp_dur_all =[];
current_lp_latency_all =[];
%need to test that opto_data and zone_data are the same length
%sometimes the zone has one extra
if length(zone_data_raw) >length(opto_data_raw)
    zone_data_raw(end,:) = [];
end

opto_data = [opto_data_raw; NaN NaN];
zone_data = [zone_data_raw ; NaN NaN];

opto_data(1,:)=[];
zone_data(1,:)=[];

opto_entry_exit = [opto_data zone_data(:,2)];

%Use nearest point to find the
%nearest exit for each opto entrance
[nearestidx nearestdistance] = nearestpoint(opto_entry_exit(:,1),opto_entry_exit(:,3),'next');

opto_entry_exit_cleaned = [opto_entry_exit(1:end-1,1) opto_entry_exit(nearestidx(1:end-1),3)];
opto_entry_exit_cleaned = [opto_entry_exit_cleaned opto_entry_exit(1:end-1,2)];
opto_entry_exit_cleaned = rmmissing(opto_entry_exit_cleaned);

% Subtract logicals from each other in Lever Press column
switches_LP_raw = diff(beh_data_raw(:,3));

LP_ON_index_raw = find(switches_LP_raw == 1)+1;
LP_ON_timestamps_raw = beh_data_raw(LP_ON_index, 1);
LP_ON_timestamps_raw = [LP_ON_timestamps_raw; NaN];

%% create an indicator just for zone entrance
%and to only look at stim within zone
zone_only_stim_indicator = NaN(length(HD_time),1);
zone_indicator = zeros(length(HD_time),1);
zone_indicator_first_vs_others =  NaN(length(HD_time),1);

for q = 1:length(opto_entry_exit_cleaned)-2
    z=1;
    while z < length(LP_ON_timestamps_raw)
        
    if LP_ON_timestamps_raw(z) >= opto_entry_exit_cleaned(q,1) && LP_ON_timestamps_raw(z) <= opto_entry_exit_cleaned(q,2)
        
        zone_indicator(z) = 1;
        
        if opto_entry_exit_cleaned(q,3) == 0
      zone_only_stim_indicator(z) = 1;
        end
        
    if opto_entry_exit_cleaned(q,3) == 1
          zone_only_stim_indicator(z) = 2;
      end
    
    
      if LP_ON_timestamps_raw(z+1) >= opto_entry_exit_cleaned(q,1) && LP_ON_timestamps_raw(z+1) <= opto_entry_exit_cleaned(q,2)
         
       z = z + find(LP_ON_timestamps_raw(z+1:end) > opto_entry_exit_cleaned(q,2),1);
               
       if z > length(LP_ON_timestamps_raw) 
           z= length(LP_ON_timestamps_raw)-1;
       end
       
      end
       
    end
     %need to skip to the next opto/entry iteration if there is a press in
    %between entry/exit. Otherwise ALL presses made within the zone would
    %betagged with a 1 or 2
    z=z+1;
    end  
end
        
n1_zone_only_stim_indicator =[NaN; zone_only_stim_indicator(1:end-1)]; 
n2_zone_only_stim_indicator =[NaN; n1_zone_only_stim_indicator(1:end-1)]; 
n3_zone_only_stim_indicator =[NaN; n2_zone_only_stim_indicator(1:end-1)]; 
n4_zone_only_stim_indicator =[NaN; n3_zone_only_stim_indicator(1:end-1)]; 
n5_zone_only_stim_indicator =[NaN; n4_zone_only_stim_indicator(1:end-1)]; 
n6_zone_only_stim_indicator =[NaN; n5_zone_only_stim_indicator(1:end-1)]; 

n_back_stim_zone_only = [zone_only_stim_indicator n1_zone_only_stim_indicator n2_zone_only_stim_indicator n3_zone_only_stim_indicator n4_zone_only_stim_indicator n5_zone_only_stim_indicator n6_zone_only_stim_indicator];

%% LME
%want to build LMEs from this data. To do this, we want to track lps and
%whether stimulation occurred prior to a given lever press, counting only
%the first lp that this happens for. Only track lps and stimulation for
%now, add in other stuff if there is anything interesting here.
%Possibl ways to do this:
%1- put an arbitrarry time limit, like if a lp occurs witihn 5s of
%stimulation we count it as stimulated
%2- don't put any time requirement, just the first lp after stimulation tag it
%3- similar to 2, but make sure that if another stimulation occurs first

%get index of the nearest lp to any given stim
closestlp = nearestpoint(Stim_ON_index, LP_ON_index,'next');

%then tag all those lps with a 1, all others with a zero
LP_stim_indicator= LP_ON_index;
closestlp =rmmissing(closestlp);
LP_stim_indicator(closestlp) = 1;
LP_stim_indicator(setdiff(1:end,closestlp)) = 0;

%we have n duration and n stimulation, now we want to get n -1 : n -10
%duration and stimulation. Could just pad with nans. 

n_back_Lengths = [];
n_lengths_array = [];
Logical_n_lengths_array= [];
Logical_n_back_Lengths =[];
press_indices = 1:length(HD_time);
press_indices = press_indices';      
n_back_stims = [];
n_stims_array = [];
n_back_stims_and_zones =[];
n_stims_and_zones_array =[];
           %for n-backs 1 to 10
       for x = 1:10 
   
       for n_it = 1:length(HD_time)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           n_back_stim = NaN;
           n_back_stim_and_zone = NaN;
           elseif  press_indices(n_it) > x
           n_back = HD_time(n_it- x);    
         n_back_stim = LP_stim_indicator(n_it - x);
         n_back_stim_and_zone = stim_and_zone_indicator(n_it - x);
           end
           n_lengths_array = [n_lengths_array n_back];
            n_stims_array = [n_stims_array n_back_stim];
            n_stims_and_zones_array = [n_stims_and_zones_array n_back_stim_and_zone];
%          %add in logical lengths 0/1 for failure/success
           Logical_n_lengths_array = n_lengths_array >= Criteria;
                     
       end
            
       n_back_Lengths = [n_back_Lengths; n_lengths_array]; 
       n_lengths_array =[];      
               
      %the logical test converts all NaNs to 0, need to put them back so we
      %don't inappropriately put lever press prior to the start of the
      %session
       Logical_n_lengths_array = double(Logical_n_lengths_array);
       Logical_n_lengths_array(1:x) = NaN;
       Logical_n_back_Lengths = [Logical_n_back_Lengths; Logical_n_lengths_array];  
       Logical_n_lengths_array =[];
       
       n_back_stims = [n_back_stims; n_stims_array]; 
       n_stims_array =[];   
       
       n_back_stims_and_zones = [n_back_stims_and_zones; n_stims_and_zones_array]; 
       n_stims_and_zones_array =[]; 
       
       end
       
 n_back_Lengths = n_back_Lengths';      
 Logical_n_back_Lengths = Logical_n_back_Lengths';     
 n_back_stims = n_back_stims'; 
 n_back_stims_and_zones =n_back_stims_and_zones'; 
 
%% Build Regression Table with behavioral variables, including stimulation indicators 
 
Regression_matrix = [HD_time n_back_Lengths Logical_n_back_Lengths stim_and_zone_indicator n_back_stims_and_zones n_back_stim_zone_only];
       
Regression_table =array2table(Regression_matrix,'VariableNames',{'n_duration','n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All', 'n_minus_four_Durations_All',...
    'n_minus_five_Durations_All', 'n_minus_six_Durations_All', 'n_minus_seven_Durations_All',...
    'n_minus_eight_Durations_All', 'n_minus_nine_Durations_All', 'n_minus_ten_Durations_All',...
    'n_minus_one_All',...
    'n_minus_two_All', 'n_minus_three_All', 'n_minus_four_All',...
    'n_minus_five_All', 'n_minus_six_All', 'n_minus_seven_All',...
    'n_minus_eight_All', 'n_minus_nine_All', 'n_minus_ten_All',...
    'n_Stim_All',...
    'n_minus_one_Stim_All','n_minus_two_Stim_All', 'n_minus_three_Stim_All', 'n_minus_four_Stim_All',...
    'n_minus_five_Stim_All', 'n_minus_six_Stim_All', 'n_minus_seven_Stim_All',...
    'n_minus_eight_Stim_All', 'n_minus_nine_Stim_All', 'n_minus_ten_Stim_All',...
    'n_stim_zone_all','n1_stim_zone_all','n2_stim_zone_all','n3_stim_zone_all','n4_stim_zone_all',...
    'n5_stim_zone_all','n6_stim_zone_all'});

%% moving average
ma_length = 60;        
moving_average_lp_length = movmean(HD_time,[ma_length,0]);
moving_average_lp_length = [NaN; NaN; moving_average_lp_length];
moving_average_lp_length = moving_average_lp_length(1:end-2);
moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; moving_average_lp_length];
moving_average_lp_length_n7andback = moving_average_lp_length_n7andback(1:end-5);
Regression_table.moving_average_lp_length_n7andback = moving_average_lp_length_n7andback;

%% IPI 
ipi1 = [NaN; diff(LP_ON_timestamps)];
ipi2 =[NaN; ipi1(1:end-1)];
Regression_table.ipi1 = ipi1;
Regression_table.ipi2 = ipi2;

%break up ipi into whether or not sitmulation occurred
ipi_and_Stim = [ipi1 LP_stim_indicator];
LP_stim_indicator = logical(LP_stim_indicator);

ipi1_stim = ipi1(LP_stim_indicator);
ipi1_nostim = ipi1(~LP_stim_indicator);
ipi1_stim_mean = nanmean(ipi1_stim);
ipi1_nostim_mean =nanmean(ipi1_nostim);
ipi1_stim_median = nanmedian(ipi1_stim);
ipi1_nostim_median =nanmedian(ipi1_nostim);

Vid_Data.ipi1_stim = ipi1_stim;
Vid_Data.ipi1_nostim = ipi1_nostim;
Vid_Data.ipi1_stim_mean = ipi1_stim_mean;
Vid_Data.ipi1_nostim_mean = ipi1_nostim_mean;
Vid_Data.ipi1_stim_median = ipi1_stim_median;
Vid_Data.ipi1_nostim_median = ipi1_nostim_median;

%% %find rewarded lps
Reward_idx = HD_time >= Criteria;
Reward_idx_n1 = [NaN; Reward_idx(1:end-1)];
% Rtotal cumulative rewards across session
reward_indicator = 0;
for q = 1:length(HD_time)
        if HD_time(q) >=Criteria
    reward_indicator = reward_indicator+1;
        end
end
 
Regression_table.Reward_idx_n1 = Reward_idx_n1; 

%% he and ipi
LP_start_times = LP_ON_timestamps;
LP_start_times =[ LP_start_times repmat(10,[length(LP_start_times) 1])  ];
 
HE_start_times = HE_ON_timestamps;
HE_start_times =[ HE_start_times repmat(12,[length(HE_start_times) 1])  ];

Event_start_times = [LP_start_times; HE_start_times ];
[spot og_idx ] = sort(Event_start_times(:,1));
sorted_Events = Event_start_times(og_idx,:);

%now loop through to find when a lp 
nanpad = [NaN NaN];

sorted_Events = [sorted_Events; nanpad];
HE_Indicator =[];
HE_counter_general =0;
for he_it = 1:length(sorted_Events)-1
              
if sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ==12  %if lp and rewarded
              
        HE_Indicator = [HE_Indicator; 1];
        
    elseif sorted_Events(he_it,2) == 10 && sorted_Events(he_it+1,2) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
end         


HE_n1_Indicator = [NaN; HE_Indicator];
HE_n1_Indicator(end) =[];
Regression_table.HE_n1_Indicator = HE_n1_Indicator;

%we have lp timestamps, but they are in bonsai time. Basically just need to
%make them in relation to the session.
%we have the start of session time (start_session_ts), so that will become time 0, and and
%then lp timestamps will be reported in relation to that
lp_on_times  = LP_ON_timestamps;
lp_on_times = lp_on_times - start_session_ts;
Regression_table.lp_on_times = lp_on_times;

%add in latency to press form stimulation 
Regression_table.stim_and_zone_latency = stim_and_zone_latency;

%% 
 %now lets create a table with n and n back lengths (actual and logical)
%and n and n back stims, and a marker for mouse and day

criteria_percent = sum(Reward_idx)/length(HD_time);
criteria_percent_indicator =  ones(length(HD_time),1)*criteria_percent;
reps = length(HD_time);
mouse_indicator = repmat(onlyID ,reps,1);
day_indicator = repmat(day_indicator,reps,1); 

%main stimulation indicator
Regression_table.stim_and_zone_indicator =stim_and_zone_indicator;
%indiactor for just first lp in zone vs others
Regression_table.zone_indicator = zone_indicator;
Regression_table.criteria_percent_indicator = criteria_percent_indicator;
Regression_table.mouse_indicator = mouse_indicator;
Regression_table.day_indicator = day_indicator;
%save the table and var names
data_table_variables = Regression_table.Properties.VariableNames;
regression_cell= table2cell(Regression_table);

Vid_Data.regression_cell = regression_cell;
Vid_Data.data_table_variables = data_table_variables;
Vid_Data.Regression_table = Regression_table;
Vid_Data.beh_data = beh_data;
Vid_Data.opto_data = opto_data;
Vid_Data.zone_data = zone_data;
