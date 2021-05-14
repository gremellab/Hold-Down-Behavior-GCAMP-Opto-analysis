function [GCAMP] = GCAMP_data_extract_r(GCAMP)
%GCAMP_data_extract_new : Saves GCAMP and Behavior data into a data structure
% beh_data : imported behavior data from 'analoginputs' .csv
% criteria : lever press length criteria in milliseconds
% mouseID : preferred naming convention for mouse ID
% training_day : preferred naming convention for training day
%% Extract behavior timestamps from analog_inputs
%in beh_data, col 2 = lever press, col 3 = HE, col 4 = reward

% Subtract logicals from each other in Lever Press column
switches_LP = diff(GCAMP.beh_data(:,2)); %2
% Subtract logicals from each other in Reinforcement column
switches_RE = diff(GCAMP.beh_data(:,4)); %4
% Subtract logicals from each other in HeadEntries column
switches_HE = diff(GCAMP.beh_data(:,3)); %3
% When the difference is 1, that's event ON
LP_ON_index = find(switches_LP == 1)+1;
RE_ON_index = find(switches_RE == 1)+1;
HE_ON_index = find(switches_HE == 1)+1;
% When the difference is -1, that's event OFF
LP_OFF_index = find(switches_LP == -1)+1;
RE_OFF_index = find(switches_RE == -1)+1;
HE_OFF_index = find(switches_HE == -1)+1;
% Get the time stamps from those event indices
LP_ON_timestamps = GCAMP.beh_data(LP_ON_index, 1);
LP_OFF_timestamps = GCAMP.beh_data(LP_OFF_index, 1);

% code to exclude the first lp or two if a lp happens too early during
% baseline (meaning there isn't enough time to have a pre onset baseline)
% LP_ON_timestamps = LP_ON_timestamps(2:end);
% LP_OFF_timestamps = LP_OFF_timestamps(2:end);

RE_ON_timestamps = GCAMP.beh_data(RE_ON_index, 1);
RE_OFF_timestamps = GCAMP.beh_data(RE_OFF_index, 1);
HE_ON_timestamps = GCAMP.beh_data(HE_ON_index, 1);
HE_OFF_timestamps = GCAMP.beh_data(HE_OFF_index,1);
% Matrix of times where Column 1 is event ON, Column 2 is event OFF
LP_All = [LP_ON_timestamps LP_OFF_timestamps];
RE_All = [RE_ON_timestamps RE_OFF_timestamps];
HE_All = [HE_ON_timestamps HE_OFF_timestamps];
% Hold down time for each lever press
HD_time=(LP_All(:,2)-LP_All(:,1))-20;
%Find Durations greater than 10s
longidx =find(HD_time>=10000);
%exclude these presses. Need to remove them from all LP variables
%(idx, start,stop ts, dur)
LP_ON_index(longidx) = [];
LP_OFF_index(longidx) = [];
LP_ON_timestamps(longidx) = [];
LP_OFF_timestamps(longidx) = [];
LP_All(longidx) = [];
HD_time(longidx) = [];

%exclude any presses that are less than 1ms after the 20ms subtraction
shortidx = find(HD_time <1);
LP_ON_index(shortidx) = [];
LP_OFF_index(shortidx) = [];
LP_ON_timestamps(shortidx) = [];
LP_OFF_timestamps(shortidx) = [];
LP_All(shortidx) = [];
HD_time(shortidx) = [];

%HE duration
HE_time =(HE_All(:,2)-HE_All(:,1))-20;

% Lever presses that met criteria and were reinforced
Criteria_met = find(HD_time >=  GCAMP.Criteria); %GCAMP.Criteria);
Criteria_fail = find(HD_time <   GCAMP.Criteria);
Total_Reinforcers = length(Criteria_met);
%HEs that were above or below the mean
HE_Above_median = find(HE_time >= median(HE_time));
HE_Below_median = find(HE_time < median(HE_time));
HE_idx = find(HE_time > 0);
% For peri-event data
SR = 20; % sampling rate

%%
%Find the first HE made after a Reinforcer
First_HE_After_RE_idx =[];
for k = 1:length(RE_ON_timestamps) %for each re_ON timestamp
    % Find index of this re_ON timestamp in the vector of all of your
    % timestamps
    Closest_HE_After_RE_idx = nearestpoint(RE_ON_timestamps(k),HE_ON_timestamps,'next');
First_HE_After_RE_idx = [First_HE_After_RE_idx; Closest_HE_After_RE_idx];
% First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
end
% Remove NaNs (HE happening outside LP)
First_HE_After_RE_idx(isnan(First_HE_After_RE_idx)) = [];
%% Save Data
% Create GCAMP data structure and save beh data
GCAMP.HoldDown_times = HD_time; %
GCAMP.Total_Reinforcers = Total_Reinforcers;
GCAMP.SR = SR;
GCAMP.Criteria_met = Criteria_met; %
GCAMP.Criteria_fail = Criteria_fail; %
GCAMP.LP_ON_timestamps = LP_ON_timestamps;
GCAMP.RE_ON_timestamps = RE_ON_timestamps;
GCAMP.LP_OFF_timestamps = LP_OFF_timestamps;
GCAMP.HE_ON_timestamps = HE_ON_timestamps;
GCAMP.HE_OFF_timestamps = HE_OFF_timestamps;
GCAMP.HE_Times = HE_time;
GCAMP.Short_HE = HE_Below_median;
GCAMP.Long_HE = HE_Above_median;
GCAMP.First_HE_After_RE =First_HE_After_RE_idx;
GCAMP.HE_idx = HE_idx;
end

