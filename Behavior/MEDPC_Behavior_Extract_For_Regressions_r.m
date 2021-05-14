function [Data] = MEDPC_Behavior_Extract_For_Regressions_BIC()
%BEHAVIOR
% The goal of this program is to take raw MEDPC data and produce usable
% information on the timing/duration of each event type (lever presses, head
% entries, reinforcements) and the intervals between events. There is
% optional code that is commented out for more variables or figures.

%When Running the code, it will prompt you to point it towards your med pc
%data. This should be formatted in 1 overarching folder per group (e.g. a
%folder containing all control or all exp. animals). Within each folder,
%should then be individual folders for each day of training. This function
%will then loop through first, the individual medpc files in the first folder
%(make sure they are named sequentially and in a standardized fashion
%across days), and then move to the next folder, where it will again loop
%through the individual medpc files, and so on.

%The end result will be saving a data structure (called Data) in the PathName_Master
%folder that contains all the relevant beh data generated from this
%function. This data file can then be used to create various models - esp.
%the linear mixed effect models in the LME_Modeling_r Script. 
%The data structure will be organized across days and mice

%This was last updated 5/10/2021 by  Drew Schreiner.
%(Gremel Lab, UCSD Dept of Psychology. Support: Please contact dcschrei at
% ucsd.edu)

tic %time the function
clear all
close all
format long %make sure matlab doesn't truncate your numbers or use scientific notation

%% -----Specify what behavior you want to extract by letter assignment with MEDPC variables as reference (i.e. LP, HE, REIN)----- %%
Alphabet = 'ABDEFGHIJKLMNOPQRSTUVWXYZC';%modified alphabet because for some reason MEDPC puts C last
LP_var = 'Q'; HE_var = 'B'; REIN_var = 'D'; Hold_Criteria_var = 'L'; %specify MEDPC variable (based on the medpc text files) for your behavior
LP_idx = strfind(Alphabet, LP_var); HE_idx = strfind(Alphabet, HE_var); REIN_idx = strfind(Alphabet, REIN_var); % index of behaviors in MEDPC file
Hold_Criteria_idx = strfind(Alphabet, Hold_Criteria_var);

%% -----Set path of your folders containing MEDPC files----- %%
Experiment_name = '33-3-beh-revised-10shuf'; % Name of experiment (will be named of saved .mat for all of the following)
PathName_Master = 'C:\Users\dcschrei\Desktop\laptop transfer\saved regression data';
PathName_Folders = uigetdir(PathName_Master); % Path with MEDPC Folder. Click on 'Select Folder' DO NOT DOUBLE CLICK
MedPCFolders = dir(PathName_Folders); % get names of MEDPC Folders
MedPCFolders = MedPCFolders(3:end); %remove invalid directory entries '.' and '..'
cd(PathName_Folders) %change directory to where MEDPC Folders are

%% -----Loop through each folder (i.e. experiment day) and within it, each subject to extract event (i.e. timing) and behavior data (i.e. # LPs)----- %%
for day = 1:length(MedPCFolders) % For each MEDPC Folder (day)
    cd([MedPCFolders(day).folder '\' MedPCFolders(day).name]) % Change directory to that folder
    PathName_Folder = cd; % get path of current directory (folder)
    MedPCFiles = dir(PathName_Folder); % get names of individual MEDPC files within that folder
    MedPCFiles = MedPCFiles(3:end); %remove invalid directory entries '.' and '..'
    for mouse = 1:length(MedPCFiles) % for each MEDPC file in this folder (mouse)
        MedPCFileName = MedPCFiles(mouse).name % get MEDPC filename
        %this is another funciton, an auto-generated textfile import from
        %matlab
        subject_vars = importfile_MEDPC_r(MedPCFileName, 11, 35); % Import first chunk of MEDPC file (which has the variables like total lever presses)
        subject_vars = subject_vars(:,1); % remove NaNs in matrix
        %--Get data from first chunk of MEDPC file variables--%
        Data.Day(day).Mouse(mouse).Name = MedPCFileName;
        Data.Day(day).Mouse(mouse).Total_Lever_Presses = subject_vars (LP_idx); % get Lever presses
        Data.Day(day).Mouse(mouse).Total_Head_Entries = subject_vars (HE_idx); % get Head entries
        Data.Day(day).Mouse(mouse).Total_Reinforcers_Earned = subject_vars (REIN_idx); % get Reinforcers earned
        Data.Day(day).Mouse(mouse).Hold_Criteria = subject_vars(Hold_Criteria_idx) * 10; % Criteria in milliseconds
        Data.Day(day).Mouse(mouse).Events_Matrix = importfile_MEDPC_r(MedPCFileName, 37, 500000); % Import event and timing data, the array in medpc files
        subject_beh = Data.Day(day).Mouse(mouse).Events_Matrix; % MEDPC Event and Timing matrix of current subject
        %% Get Timestamps of Each Event
        time = reshape(((floor(subject_beh))'),[],1); % Get the integer part of the number (ie: the number to the left of the . which is a timestamp relative to the previous event)
        time(isnan(time)) = []; % Remove NaNs
        time = cumsum(time)*10; % Cumulative sum of event timestamps to get the time in session instead of relative to the preceding event
        event = (reshape(((subject_beh-floor(subject_beh))'),[],1)) * 100; %Get the right of the decimal part of the number (ie: the evnent code, e.g. 10 = lever press start)
        event(isnan(event)) = []; % Remove NaNs
        time_with_events = cat(2,event,time);   % Concatenate time and events
        time_with_events(any(isnan(time_with_events), 2), :) = []; % Remove NaNs
        time_with_events(:,1) = round(time_with_events(:,1)); % Round all the event types to remove decimal notation
        time_with_events = unique(time_with_events,'rows','stable'); % Remove duplicates due to med-pc time resolution
        time_with_events(find(time_with_events(:,1) == 30),:) = []; % Removes event code 30
            
        %find the last time the .15 variable, which is lever press end, add 1 to include the reinforcer .20 if it is there. 
        %This excludes the session end variable, or other spurious events after the final press 
       lastleverstop =  find(time_with_events(1:end,1)==31,1,'last');
       eventlengthARRAY = []; j = 1; 
        
                while j<lastleverstop  %while loop to go through the event array up until just before the 31 session end variable
%            %loop through all the events and find how much time elapses
%            between pairs of events
             eventtime = time_with_events(j,2);
             eventafter = time_with_events(j+1,1);
             eventaftertime = time_with_events(j+1,2);
             eventlength = (eventaftertime-eventtime);
             eventlengthARRAY = [eventlengthARRAY eventlength];
            j = j +1;
               end
        
          time_with_events(find(time_with_events(:,1) == 31,1,'last')+1:end,:) = []; % Removes anything after event code 31
            %look for errors of alignment
        try
            time_with_events = [time_with_events [(eventlengthARRAY)'; NaN]];
        catch
            time_with_events = [time_with_events [(eventlengthARRAY)']];
        end
                
         %for loop to exlude any times there is a lever press that is not
         %paired with a stop variable, due to med-pc time resolution, or end of session errors          
        for z = 1:length(time_with_events)
            if time_with_events(z,1) ==10 %LP code
                if isempty(find(time_with_events(z+1) == 15 || time_with_events(z+1,1) ==35)) %if there is not a LP stop(15) or LP Stimulation code (35)
                    time_with_events(z,:) = NaN; %replace with NaN 
                    event(z,:) = NaN;
                end
            end
        end       
        %now go in the reverse direction, make sure each 15 stop is
        %preceded by either LP Start (10) or LP Stimulation (35)
               for z = 1:length(time_with_events)
            if time_with_events(z,1) ==15  %LP stop code
                if isempty(find(time_with_events(z-1) == 10  || time_with_events(z-1) ==35))
                    time_with_events(z,:) = NaN;
                    event(z,:) = NaN;
                end
            end
        end   
           
        %% This code is needed for Optogenetic Stimulation sessions.
        %the 35 stimulation variable occurs after the start but before the
        %stop. Need to take this into account to get the events array
        %timing properly aligned for lever press durations.
        %look for the 35 stim variable, and paste the time associated with it into the time
        %with events array aligned with the 10 (to get duration of LP,
        %since the 35 occurs in between 10 and 15)
        time_with_events_dummy = time_with_events;
        for lp_it = 1:length(time_with_events-1)
            if time_with_events(lp_it,1) == 10 &&  time_with_events(lp_it+1,1)==35
                time_with_events_dummy(lp_it,3) = time_with_events_dummy(lp_it+1,3);
            end
        end
        
                 %% Removal of Lever Presses over 10000ms (10s)  
           for q = 1:length(time_with_events_dummy-1)
                if time_with_events_dummy(q,1) ==10 && time_with_events_dummy(q,3) >=10000 
                    
                %turn them to Nans (can't remove or the loop will break as
                %the length will change)
                time_with_events_dummy(q,:) =[NaN];
                time_with_events_dummy(q+1,:) =[NaN];
                %special if for probabalstic reward, since a long press will
                %not always be followed by the 20 reward variable. Don't
                %want to remove the next event of it is unrewarded
                if time_with_events(q+2) == 20
                  time_with_events_dummy(q+2,:) =[NaN];
                end
                %do the same to the event and time single arrays
                time(q) =[NaN];
                time(q+1) = [NaN];
                 time(q+2) = [NaN];
                event(q)= [NaN];
                event(q+1) = [NaN];
                event(q+2) = [NaN];
                end
                
                %also exclude if press is 0s after subtracting 20ms for
                %mepc time
                if time_with_events_dummy(q,1) ==10 && time_with_events_dummy(q,3) <= 20 
                    
                %turn them to Nans (can't remove or the loop will break as
                %the length will change)
                time_with_events_dummy(q,:) =[NaN];
                time_with_events_dummy(q+1,:) =[NaN];
                %special if for probabalstic reward, since a long press will
                %not always be followed by the 20 reward variable. Don't
                %want to remove the next event of it is unrewarded
                if time_with_events(q+2) == 20
                  time_with_events_dummy(q+2,:) =[NaN];
                end
                %do the same to the event and time single arrays
                time(q) =[NaN];
                time(q+1) = [NaN];
                time(q+2) = [NaN];
                event(q)= [NaN];
                event(q+1) = [NaN];
                event(q+2) = [NaN];
                end
                
                         end
                                     
                time_with_events = time_with_events_dummy;    
      
 %% Lever press timing and durations               
        Lever_Press_Matrix = time_with_events(find(time_with_events(:,1) == 10),:); % Matrix of event code, timestamp, and length
        Lever_Press_idx = find(time_with_events(:,1) == 10); % Index of when Lever Press Ocurred
        Lever_Press_Lengths = Lever_Press_Matrix(:,3)-20; % Subtract 20 to get true lever press lengths (due to medpc timing)             
        Lever_Press_True_End = Lever_Press_Matrix(:,2) + Lever_Press_Matrix(:,3); % Lever Press Stop Timestamp
        Lever_Press_IPI = Lever_Press_Matrix(2:end,2) - Lever_Press_True_End(1:end-1); % Lever Press Inter Press Interval
        Lever_Press_ts = Lever_Press_Matrix(:,2); % Timestamps of when Lever Press Ocurred in total session time
        %add a 0 to pad the IPI array so it is the same length as lever
        %presses (since there can be no IPI for the first press)
        Lever_Press_IPI = [NaN;Lever_Press_IPI];
        Data.Day(day).Mouse(mouse).Lever_Press_ts = Lever_Press_ts;
        Data.Day(day).Mouse(mouse).Lever_Press_duration = Lever_Press_Lengths;
        Data.Day(day).Mouse(mouse).Lever_Press_IPI = Lever_Press_IPI;
        Data.Day(day).Mouse(mouse).Lever_Press_idx = Lever_Press_idx;
        Data.Day(day).Mouse(mouse).Events_Time_Matrix = time_with_events; % Add to Data structure 
        
%% Lever Press Performance
        Lever_Press_Rewarded_idx = find(Lever_Press_Lengths >= Data.Day(day).Mouse(mouse).Hold_Criteria); % Index of when successful Lever Press Ocurred
        Lever_Press_Rewarded_Lengths = Lever_Press_Lengths(Lever_Press_Rewarded_idx); % Length successful Lever Press
        Lever_Press_Rewarded_ts = Lever_Press_Matrix(Lever_Press_Rewarded_idx,2);
        
        %check for if there were no rewarded presses
        if isempty(Lever_Press_Rewarded_idx)
            Lever_Press_Rewarded_idx = NaN;
            Lever_Press_Rewarded_Lengths = NaN;
            Lever_Press_Rewarded_ts = NaN;
        end
        
        Lever_Press_Not_Rewarded_idx = find(~(Lever_Press_Lengths >= Data.Day(day).Mouse(mouse).Hold_Criteria)); % Index of when unsuccessful Lever Press Ocurred
        Lever_Press_Not_Rewarded_Lengths = Lever_Press_Lengths(Lever_Press_Not_Rewarded_idx); % Length unsuccessful Lever Press
        Lever_Press_Not_Rewarded_ts = Lever_Press_Matrix(Lever_Press_Not_Rewarded_idx,2);
        
        %check for if there were no unrewarded presses
         if isempty(Lever_Press_Not_Rewarded_idx)
            Lever_Press_Not_Rewarded_idx = NaN;
            Lever_Press_Not_Rewarded_Lengths = NaN;
            Lever_Press_Not_Rewarded_ts = NaN;
         end
        
        %save in structure
        Data.Day(day).Mouse(mouse).Lever_Press_Rewarded_idx = Lever_Press_Rewarded_idx;
        Data.Day(day).Mouse(mouse).Lever_Press_Rewarded_Lengths = Lever_Press_Rewarded_Lengths;
        Data.Day(day).Mouse(mouse).Lever_Press_Rewarded_ts = Lever_Press_Rewarded_ts;
        Data.Day(day).Mouse(mouse).Lever_Press_Not_Rewarded_idx = Lever_Press_Not_Rewarded_idx;
        Data.Day(day).Mouse(mouse).Lever_Press_Not_Rewarded_Lengths = Lever_Press_Not_Rewarded_Lengths;
        Data.Day(day).Mouse(mouse).Lever_Press_Not_Rewarded_ts = Lever_Press_Not_Rewarded_ts;

        %% Feedback, i.e., change in duration after a rewarded or unrewarded lever press
        Lever_Press_Length_After_Success_idx = Lever_Press_Rewarded_idx + 1; % Look for lever press index after success
        if isnan(Lever_Press_Length_After_Success_idx) %check to see if there is a press after success, otherwise NaN
            Lever_Press_Length_After_Success_remove_idx = NaN;
        else
        Lever_Press_Length_After_Success_remove_idx = find(Lever_Press_Length_After_Success_idx > length(Lever_Press_ts)); % Remove indices outside total lever presses
        end
              
        Lever_Press_Length_After_Failure_idx = Lever_Press_Not_Rewarded_idx + 1; % Look for lever press index after failure
          if isnan(Lever_Press_Length_After_Failure_idx) %check to see if there is a press after failure otherwise NaN
              Lever_Press_Length_After_Failure_remove_idx = NaN;
          else
        Lever_Press_Length_After_Failure_remove_idx = find(Lever_Press_Length_After_Failure_idx > length(Lever_Press_ts)); % Remove indices outside total lever presses
          end
                 
        Lever_Press_Length_After_Success = Lever_Press_Lengths(Lever_Press_Length_After_Success_idx(1:end-length(Lever_Press_Length_After_Success_remove_idx))); % Lever Press Length at n+1 after success
        Lever_Press_Length_After_Failure = Lever_Press_Lengths(Lever_Press_Length_After_Failure_idx(1:end-length(Lever_Press_Length_After_Failure_remove_idx))); % Lever Press Length at n+1 after failure
        Lever_Press_Length_After_Success_Change = Lever_Press_Length_After_Success... % Lever Press Duration change at n+1 after success
            - Lever_Press_Lengths(Lever_Press_Rewarded_idx(1:end-length(Lever_Press_Length_After_Success_remove_idx)));
        Lever_Press_Length_After_Failure_Change = Lever_Press_Length_After_Failure...% Lever Press Duration change at n+1 after failure
            - Lever_Press_Lengths(Lever_Press_Not_Rewarded_idx(1:end-length(Lever_Press_Length_After_Failure_remove_idx)));
       %check for press after failure or success being emppty
       
        if isempty(Lever_Press_Length_After_Success)
            Lever_Press_Length_After_Success = NaN;
            Lever_Press_Length_After_Success_Change = NaN;
        end
        if isempty(Lever_Press_Length_After_Failure)
            Lever_Press_Length_After_Failure = NaN;
            Lever_Press_Length_After_Failure_Change = NaN;
        end
        
        %save in structure
        Data.Day(day).Mouse(mouse).Lever_Press_Length_After_Success = Lever_Press_Length_After_Success;
        Data.Day(day).Mouse(mouse).Lever_Press_Length_After_Failure = Lever_Press_Length_After_Failure;
        Data.Day(day).Mouse(mouse).Lever_Press_Length_After_Success_Change = Lever_Press_Length_After_Success_Change;
        Data.Day(day).Mouse(mouse).Lever_Press_Length_After_Failure_Change = Lever_Press_Length_After_Failure_Change;
           
%% Moving mean of LP Length
       % want to calcuulate the moving mean of durations, NOT include the current lp or n-1 
       ma_length = 60; %length of ma average, arrived at using BIC for different lengths
       moving_average_lp_length = movmean(Lever_Press_Lengths,[ma_length,0]);
       %simply pad the start with 2 Nans, this will then give you the 60th
       %back moving average from n-2>n-60 aligned with the appropriate n
       moving_average_lp_length = [NaN; NaN; moving_average_lp_length];
       moving_average_lp_length = moving_average_lp_length(1:end-2);
       Data.Day(day).Mouse(mouse).moving_average_lp_length = moving_average_lp_length;
       
       %Instead, exclude presses n - 1 through n - 6, since those are individual related to press n
       moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; moving_average_lp_length];
       moving_average_lp_length_n7andback = moving_average_lp_length_n7andback(1:end-5);
       Data.Day(day).Mouse(mouse).moving_average_lp_length_n7andback = moving_average_lp_length_n7andback;
 
%% n-x lever press lengths
%create a matrix where we can align press n with presses n - x in separate
%columns
%create variables
           n_back_Lengths = {};
           n_lengths_array = [];
           Logical_n_lengths_array= [];
           Logical_n_back_Lengths ={};
           press_indices = 1:length(Lever_Press_Lengths);
           press_indices = press_indices';
           n_to_n_back_ipi_array =[];
           n_to_n_back_ipi_all ={};
                             
           %for looping across as many n-backs as you want (10 here)
       for x = 1:10 
   %then loop within each n-back, to get the n - x of each lever press
       for n_it = 1:length(Lever_Press_Lengths)
           %need to check that if there are at least x presses before the
           %current state press.
           if  press_indices(n_it) <= x
           %if there are not, then put a NaN there
           n_back = NaN;
           n_to_n_back_ipi = NaN;
           elseif  press_indices(n_it) > x
           n_back = Lever_Press_Lengths(n_it- x);    
           n_to_n_back_ipi = sum(Lever_Press_IPI(n_it-x:n_it)); %also get n-back ipi here
           end
           n_lengths_array = [n_lengths_array n_back];
           n_to_n_back_ipi_array = [n_to_n_back_ipi_array n_to_n_back_ipi];       
           %add in logical lengths 0/1 for failure/success
           Logical_n_lengths_array = n_lengths_array >= Data.Day(day).Mouse(mouse).Hold_Criteria;                       
       end       
      
       n_back_Lengths = [n_back_Lengths; n_lengths_array]; 
       n_lengths_array =[];      
       Logical_n_back_Lengths = [Logical_n_back_Lengths; Logical_n_lengths_array];  
       Logical_n_lengths_array =[];     
       n_to_n_back_ipi_all =[n_to_n_back_ipi_all; n_to_n_back_ipi_array ];
       n_to_n_back_ipi_array=[];
          
       end
       
       n_to_n_back_ipi_all = [Lever_Press_IPI';n_to_n_back_ipi_all];       
             
       %the logical tests turn NaNs into 0, but we need the NaNs to pad the
       %n-back array
       %loop through the array and convert the first x logical values to
       %Nans to make sure they line up with the absolute values
           for i = 1:length(Logical_n_back_Lengths)
           Logical_n_back_Lengths{i,1} =double(Logical_n_back_Lengths{i,1});
           Logical_n_back_Lengths{i,1}(1,1:i) = NaN;
           end
        
      Data.Day(day).Mouse(mouse).n_back_Lengths = n_back_Lengths;
      Data.Day(day).Mouse(mouse).Logical_n_back_Lengths = Logical_n_back_Lengths;
      Data.Day(day).Mouse(mouse).n_back_IPIs = n_to_n_back_ipi_all;
      
%% LP and Headentry (HE) Counting
%look at patterning of LPs and HEs
%loop through the event array
%find proportion of LPs (whether rewarded or unrewarded) were followed by
%HE or not - not taking time into account
event =round(event);
 Total_LPs = sum(event==15) ;%total # lps
 Total_HEs = sum(event==12); %total HEs
 Total_RE = sum(event == 20);
 Total_Unrewarded = Total_LPs-Total_RE;
 HE_counter_rewarded = 0;
 HE_counter_general = 0;
 HE_counter_unrewarded = 0;
 
 HE_Indicator = []; %will be an array to indicate if LP was followed by HE
 
 session_end = find(event==31);
 if isempty(session_end)
     session_end = length(event)-1;
 end
 %loop through event array to find when lps are followed by HE, and if this
 %is different for rewarded/unrewarded presses
for i =1:session_end
if event(i) == 15 && event(i+1) ==20  %if lp and rewarded
    
    if event(i+2) == 12 %if theres a 12 after the 20 (ie, a HE) add it to the general he counter, reward specific, 
        %and the logical HE_indicator array
        HE_counter_rewarded = HE_counter_rewarded +1;
        HE_counter_general = HE_counter_general +1;
        HE_Indicator = [HE_Indicator; 1];
        
    elseif event(i+2) ~= 12 %otherwise don't
             HE_Indicator = [HE_Indicator; 0];
    end
end

if event(i) == 15 && event(i+1) ~=20  %if LP was unrewarded
    if event(i+1) == 12 %if theres a 12 after the 15 add it to the general he counter and unreward specific
        HE_counter_unrewarded = HE_counter_unrewarded +1;
        HE_counter_general = HE_counter_general +1;
        HE_Indicator = [HE_Indicator; 1];
        
    elseif event(i+1) ~= 12
             HE_Indicator = [HE_Indicator; 0];
    end
             
end
end
%proportion of lp types followed by HE
Total_HE_Prop = HE_counter_general/Total_LPs;
Reward_HE_Prop = HE_counter_rewarded/Total_RE;
Unreward_HE_Prop = HE_counter_unrewarded/Total_Unrewarded;

 Data.Day(day).Mouse(mouse).Total_HE_Prop = Total_HE_Prop;
 Data.Day(day).Mouse(mouse).Reward_HE_Prop = Reward_HE_Prop;
 Data.Day(day).Mouse(mouse).Unreward_HE_Prop = Unreward_HE_Prop;
 Data.Day(day).Mouse(mouse).HE_Indicator = HE_Indicator; 

%add a NaN to the start of HE_indicator to make it be an indicator for the
%n-1 lp
HE_n_1_Indicator = [NaN; HE_Indicator];
HE_n_1_Indicator(end) =[];
Data.Day(day).Mouse(mouse).HE_n_1_Indicator = HE_n_1_Indicator;
%ditto for n -2
HE_n_2_Indicator = [NaN; NaN; HE_Indicator];
HE_n_2_Indicator(end-1:end) =[];
Data.Day(day).Mouse(mouse).HE_n_2_Indicator = HE_n_2_Indicator;
%n-3
HE_n_3_Indicator = [NaN; NaN; NaN; HE_Indicator];
HE_n_3_Indicator(end-2:end) =[];
Data.Day(day).Mouse(mouse).HE_n_3_Indicator = HE_n_3_Indicator;
%n-4
HE_n_4_Indicator = [NaN; NaN; NaN; NaN; HE_Indicator];
HE_n_4_Indicator(end-3:end) =[];
Data.Day(day).Mouse(mouse).HE_n_4_Indicator = HE_n_4_Indicator;
%n-5
HE_n_5_Indicator = [NaN; NaN; NaN; NaN; NaN; HE_Indicator];
HE_n_5_Indicator(end-4:end) =[];
Data.Day(day).Mouse(mouse).HE_n_5_Indicator = HE_n_5_Indicator;
%n-6
HE_n_6_Indicator = [NaN; NaN; NaN; NaN; NaN; NaN; HE_Indicator];
HE_n_6_Indicator(end-5:end) =[];
Data.Day(day).Mouse(mouse).HE_n_6_Indicator = HE_n_6_Indicator;

%add in indicator variable for overall criteria percent on a given day that's as long as their are
%lever presses. This will be used as a covariate in the model for overall
%performance
%an array that's total lever presses made length
%need to make sure that the criteria percent correctly rflects successful
%presses and not just rewards when there is probabalistic reward
Criteria_Percent = 100*(sum(Lever_Press_Lengths >=Data.Day(day).Mouse(mouse).Hold_Criteria) / length(Lever_Press_Lengths));
Data.Day(day).Mouse(mouse).Criteria_Percent = Criteria_Percent;
crit_percent_indicator = ones(length(Lever_Press_Lengths),1)*Data.Day(day).Mouse(mouse).Criteria_Percent;
criteria_indicator = ones(length(Lever_Press_Lengths),1)*Data.Day(day).Mouse(mouse).Hold_Criteria;
avg_duration_indicator = ones(length(Lever_Press_Lengths),1)*mean(Lever_Press_Lengths);
Data.Day(day).Mouse(mouse).crit_percent_indicator = crit_percent_indicator;
Data.Day(day).Mouse(mouse).criteria_indicator = criteria_indicator; %an indicator of what the criteria itself was (800ms, 1600ms)
Data.Day(day).Mouse(mouse).avg_duration_indicator = avg_duration_indicator; %avg press duration in a day

%% Patterning of Rewards
 RE_Indicator = []; 
 %add in a logical index for whether or not a press is rewarded based on
 %exceeding criteria
 time_with_events_with_reward = time_with_events(:,3) >=  Data.Day(day).Mouse(mouse).Hold_Criteria;
 time_with_events_with_reward = [time_with_events time_with_events_with_reward];
 %pad with a NAN
 time_with_events_with_reward = [time_with_events_with_reward; NaN NaN NaN NaN];
 %remove HE events (12 and 13) for ease
 rem_idx = time_with_events_with_reward(:,1) == 12 ;
 time_with_events_with_reward(rem_idx,:) =[];
  rem_idx = time_with_events_with_reward(:,1) == 13 ;
 time_with_events_with_reward(rem_idx,:) =[];
 
 session_end = find(time_with_events_with_reward(:,1)==31);
 if isempty(session_end)
     session_end = length(time_with_events_with_reward)-1;
 end
 
 %% This is for probabalistic reward. We check to see if presses that exceeded 
 %the criterion were actually rewarded (followed by a 20 reward variable)
 %then we use "2" as an indicator for over criteria + reward, a "1" for
 %over criteria - reward and a "0" for under criteria. Absent probabilstic
 %reward, it is just 0/1.
 %won't affect anything when reward is deterministic
for i =1:session_end
if time_with_events_with_reward(i,1) == 10 && time_with_events_with_reward(i,4) == 1  %if lp and rewarded
    
    if time_with_events_with_reward(i+2,1) == 20 %if theres a 20 after the lp stop variable (10>15>20), that lp was actually rewarded
        RE_Indicator = [RE_Indicator; 2];
        
    elseif time_with_events_with_reward(i+2,1) ~= 20       
             RE_Indicator = [RE_Indicator; 1];
    end
end
if time_with_events_with_reward(i,1) == 10 && time_with_events_with_reward(i,4) == 0 %also put a 0 if they were too short
    
RE_Indicator = [RE_Indicator; 0];
end

end

Data.Day(day).Mouse(mouse).RE_Indicator =RE_Indicator;
%this will tell us if n was rewarded. For n-back, simply pad with nans 
RE_n_1_Indicator = [NaN; RE_Indicator];
RE_n_1_Indicator(end) =[];
Data.Day(day).Mouse(mouse).RE_n_1_Indicator = RE_n_1_Indicator;
%n-2
RE_n_2_Indicator = [NaN; NaN; RE_Indicator];
RE_n_2_Indicator(end-1:end) =[];
Data.Day(day).Mouse(mouse).RE_n_2_Indicator = RE_n_2_Indicator;
%n-3
RE_n_3_Indicator = [NaN; NaN; NaN; RE_Indicator];
RE_n_3_Indicator(end-2:end) =[];
Data.Day(day).Mouse(mouse).RE_n_3_Indicator = RE_n_3_Indicator;
%n-4
RE_n_4_Indicator = [NaN; NaN; NaN; NaN; RE_Indicator];
RE_n_4_Indicator(end-3:end) =[];
Data.Day(day).Mouse(mouse).RE_n_4_Indicator = RE_n_4_Indicator;

%% Stim indicator
%for use when stimulation occurs optogenetically. Will not affect normal
%days.
stim_ind_all =[]; %logical indicator of if a LP is stimulated
for stim_it = 1:length(time_with_events-1);
    if time_with_events(stim_it,1) == 10 && time_with_events(stim_it+1,1) == 35 %35 is the code for stimulation. If 10 is followed by 35 = stim.
        stim_ind = 1;
        stim_ind_all = [stim_ind_all; stim_ind];
    elseif time_with_events(stim_it,1) == 10 && time_with_events(stim_it+1,1) ~= 35 
        stim_ind = 0;
        stim_ind_all = [stim_ind_all; stim_ind];
    end
end
%pad the stim indicator with Nans to get the lags
stim_ind_all_n1 = [NaN; stim_ind_all(1:end-1)];
stim_ind_all_n2 = [NaN; stim_ind_all_n1(1:end-1)];
stim_ind_all_n3 = [NaN; stim_ind_all_n2(1:end-1)];
stim_ind_all_n4 = [NaN; stim_ind_all_n3(1:end-1)];
stim_ind_all_n5 = [NaN; stim_ind_all_n4(1:end-1)];
stim_ind_all_n6 = [NaN; stim_ind_all_n5(1:end-1)];
Data.Day(day).Mouse(mouse).stim_ind_all = stim_ind_all;
Data.Day(day).Mouse(mouse).stim_ind_all_n1 = stim_ind_all_n1;
Data.Day(day).Mouse(mouse).stim_ind_all_n2 = stim_ind_all_n2;
Data.Day(day).Mouse(mouse).stim_ind_all_n3 = stim_ind_all_n3;
Data.Day(day).Mouse(mouse).stim_ind_all_n4 = stim_ind_all_n4;
Data.Day(day).Mouse(mouse).stim_ind_all_n5 = stim_ind_all_n5;
Data.Day(day).Mouse(mouse).stim_ind_all_n6 = stim_ind_all_n6;

%% State Performance
% Determine States via cumulative sum of lever press durations
%this is using the overall mean and std of the entire days data to
%calculate a rolling upper/lowersum of the durations, looking for times
%when the uppersum exceedes 2 SE relative to the mean

%need to chehck that there are any lps in order to get mean and std for
%cusum
if isempty(Lever_Press_Lengths)
    Up_State_idx = NaN;
    ilower = NaN;
    uppersum = NaN;
    lowersum = NaN;
    uppersum_sem = NaN;
    Up_State_Lengths = NaN;
    Up_State_ts = NaN;
else

    if std(Lever_Press_Lengths) == 0 %if std is 0 because of only 1 LP, just code it as 1ms to avoid error
        [Up_State_idx,ilower,uppersum,lowersum] = cusum(Lever_Press_Lengths,2,1,mean(Lever_Press_Lengths),1,'all');
        uppersum_sem = uppersum/mean(Lever_Press_Lengths);
        Up_State_Lengths = Lever_Press_Lengths(Up_State_idx); % Lever Press Lengths in upstate
        Up_State_ts = Lever_Press_ts(Up_State_idx); % Timestamp of Lever Press in upstate      
    else
        [Up_State_idx,ilower,uppersum,lowersum] = cusum(Lever_Press_Lengths,2,1,mean(Lever_Press_Lengths),std(Lever_Press_Lengths),'all');
        uppersum_sem = uppersum/mean(Lever_Press_Lengths);
        Up_State_Lengths = Lever_Press_Lengths(Up_State_idx); % Lever Press Lengths in upstate
        Up_State_ts = Lever_Press_ts(Up_State_idx); % Timestamp of Lever Press in upstate       
    end   
end

        %need to get a logical  column of 1 and 0 to indicate up-state or not
        idxxx=[];
        for state_it = 1:length(Up_State_idx)
          idxxx=[idxxx  press_indices==Up_State_idx(state_it)];
        end
        Up_State_idx_Logical=any(idxxx,2);
        %if very few presses are made, then sometimes no presses will be in
        %upstate. If so, create matrix of 0 that is as long as lps
                if isempty(Up_State_idx_Logical)
                    Up_State_idx_Logical = zeros(length(Lever_Press_Lengths),1);
                end
            
         %get the proportion of up to down state
         up_state_prop = sum(Up_State_idx_Logical)/length(Up_State_idx_Logical);
         Data.Day(day).Mouse(mouse).up_state_prop = up_state_prop;
       
%loop through the upstate index to find the avg consectuvie number of up state presses
sequence_above = 0;
conseq =0;
sequence_counter =[0];
 Up_State_idx_pad = [Up_State_idx; NaN]; %pad with one nan
for upper_it = 1:length(Up_State_idx)-1
   current = Up_State_idx(upper_it);
   if Up_State_idx(upper_it+1) == current +1
       conseq = conseq + 1;
     if conseq == 1
        sequence_above = sequence_above + 2;
       else 
        sequence_above = sequence_above + 1;
         continue
     end
       
   else
          sequence_counter = [sequence_counter conseq+1];
          conseq = 0;
       continue
  end
         
end
%add last idx on since we cant loop past the length
%remove the first 0, since it reflects the nan pad from the loop
 sequence_counter = [sequence_counter conseq+1];
 sequence_counter(sequence_counter==0) =[];
 
sequence_counter(1) = sequence_counter(1) +1;
sequence_length_mean_including_ones = mean(sequence_counter);
Data.Day(day).Mouse(mouse).Up_State_Consecutive = sequence_counter;
Data.Day(day).Mouse(mouse).Up_State_Consecutive_mean = sequence_length_mean_including_ones;
                
        %get the state identity of n-1, to be used to see if the state
        %identy of n-1 interacts with n:n-1
        Up_State_idx_n1= [NaN; Up_State_idx_Logical];
        Up_State_idx_n1=Up_State_idx_n1(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n1 = Up_State_idx_n1;
        %n-2
        Up_State_idx_n2= [NaN; Up_State_idx_n1];
        Up_State_idx_n2=Up_State_idx_n2(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n2 = Up_State_idx_n2;
        %n-3
        Up_State_idx_n3= [NaN; Up_State_idx_n2];
        Up_State_idx_n3=Up_State_idx_n3(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n3 = Up_State_idx_n3;
        %n-4
        Up_State_idx_n4= [NaN; Up_State_idx_n3];
        Up_State_idx_n4=Up_State_idx_n4(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n4 = Up_State_idx_n4;
        %n-5
        Up_State_idx_n5= [NaN; Up_State_idx_n4];
        Up_State_idx_n5=Up_State_idx_n5(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n5 = Up_State_idx_n5;
        %n-6
        Up_State_idx_n6= [NaN; Up_State_idx_n5];
        Up_State_idx_n6=Up_State_idx_n6(1:end-1);
        Data.Day(day).Mouse(mouse).Up_State_idx_n6 = Up_State_idx_n6;
        %is the up-state press also rewarded?
        Up_State_Reward_idx = intersect(Lever_Press_Rewarded_idx,Up_State_idx); % Index of up-state rewarded presses
        Up_State_Reward_Lengths = Lever_Press_Lengths(Up_State_Reward_idx); % Length of up-state rewarded presses
        Up_State_Reward_ts = Lever_Press_ts(Up_State_Reward_idx); % Timestamp of up-state rewarded presses
        %is the up-state press unrewarded?
        Up_State_Fail_idx = intersect(Lever_Press_Not_Rewarded_idx,Up_State_idx); % Index of up-state unrewarded presses
        Up_State_Fail_Lengths = Lever_Press_Lengths(Up_State_Fail_idx); % Length of up-state unrewarded presses
        Up_State_Fail_ts = Lever_Press_ts(Up_State_Fail_idx); % Timestamp of up-state unrewarded presses
        
          Data.Day(day).Mouse(mouse).Up_State_idx_Logical = Up_State_idx_Logical;
          Data.Day(day).Mouse(mouse).uppersum_sem = uppersum_sem;

%% total cumulative rewards across session
%used to look for across session/satiation effects
total_reward_indicator = [];
reward_indicator = 0;
for q = 1:length(Lever_Press_Lengths)
    total_reward_indicator = [total_reward_indicator; reward_indicator];
        if Lever_Press_Lengths(q) >= Data.Day(day).Mouse(mouse).Hold_Criteria
    reward_indicator = reward_indicator+1;
        end
end
  
%Alternative for probabalistic reward, when not every press over criteria
%is rewarded
total_reward_indicator_prob =[];
reward_indicator_prob = 0;
for k = 1:length(RE_Indicator)
        total_reward_indicator_prob = [total_reward_indicator_prob; reward_indicator_prob];
    if RE_Indicator(k) ==2
        reward_indicator_prob = reward_indicator_prob+1;
    end
end
        
    Data.Day(day).Mouse(mouse).total_reward_indicator = total_reward_indicator;
    Data.Day(day).Mouse(mouse).total_reward_indicator_prob = total_reward_indicator_prob;

%% Code to shuffle the order shuf_it number of times in order to build a distribution of shuffled data.
%also need to build the n-back arrays from these shuffled datasets
%we also want to preserve information like ipi, HE indicator, and timestamp
%for an individual press, just shuffling the order
%create a matrix that includes all these variables spread out across
%columns, then shuffle the rows only

ipi_nback_matrix =[];
for ipi_its = 1:length(n_to_n_back_ipi_all)
 ipi_nback_matrix = [ipi_nback_matrix n_to_n_back_ipi_all{ipi_its}'];
end
%create shuffled versions of all the variables we care about for LMEs
Shuffled_Durs_Distribtuion = [];
shuffled_ipi1_Distribution = [];
shuffled_ipi2_Distribution = [];
shuffled_ipi3_Distribution = [];
shuffled_he1_Distribution =[];
shuffled_he2_Distribution =[];
shuffled_he3_Distribution =[];
shuffled_he4_Distribution =[];
shuffled_MA_Distribution =[];
Shuffled_ts_Distribution =[];
Shuffled_Up_State_idx_n1_Distribution =[];
shuffled_MA_n7_Distribution =[];
shuffled_total_reward_indicator_dist =[];

    Shuffled_stim_ind_all_Distribution =[]; 
    Shuffled_stim_ind_all_n1_Distribution =[];    
    Shuffled_stim_ind_all_n2_Distribution =[]; 
    Shuffled_stim_ind_all_n3_Distribution =[];
    Shuffled_stim_ind_all_n4_Distribution =[];
    Shuffled_stim_ind_all_n5_Distribution =[];
    Shuffled_stim_ind_all_n6_Distribution =[];

mean_Shuffled_Lever_Press_Length_After_Success_Dist = [];
mean_Shuffled_Lever_Press_Length_After_Failure_Dist = [];         
mean_Shuffled_Lever_Press_Change_After_Success_Dist = [];
mean_Shuffled_Lever_Press_Change_After_Failure_Dist = [];         
std_Shuffled_Lever_Press_Length_After_Success_Dist = [];
std_Shuffled_Lever_Press_Length_After_Failure_Dist = [];   
entropy_Shuffled_Lever_Press_Length_After_Success_Dist = [];
entropy_Shuffled_Lever_Press_Length_After_Failure_Dist = [];
Up_State_Consecutive_Shuf_Dist = [];
Up_State_Consecutive_Shuf_Mean = [];
Up_State_Prop_Above_Shuf =[];    
  %shuf_it controls how many shuffled versions we create. 100 is sufficient for a first pass
  %But need to have 1000 to be more sure 
for shuf_it = 1:10
    %first we create a shuffled version of the variable by copying it over
    Shuffled_Durs = Lever_Press_Lengths;
    Shuffled_Durs = Shuffled_Durs(randperm(length(Shuffled_Durs))); %then use randperm to shuffle the order of that variable
    Shuffled_Durs_Distribtuion = [Shuffled_Durs_Distribtuion Shuffled_Durs];  %then add it to an array which will eventually have all the shuffled versions in it
    %time stamps
    Shuffled_ts = Lever_Press_ts;
    Shuffled_ts = Shuffled_ts(randperm(length(Shuffled_ts)));
    Shuffled_ts_Distribution = [Shuffled_ts_Distribution Shuffled_ts];
        
    %calculate a moving average for each shuffle
    shuffled_moving_average_lp_length = movmean(Shuffled_Durs,[ma_length,0]);      
    shuffled_moving_average_lp_length = [NaN; NaN; shuffled_moving_average_lp_length];  
    shuffled_moving_average_lp_length = shuffled_moving_average_lp_length(1:end-2);
    shuffled_MA_Distribution = [shuffled_MA_Distribution shuffled_moving_average_lp_length];
     %n7: ma_length mov. average
    moving_average_lp_length_n7andback = [NaN; NaN; NaN; NaN; NaN; shuffled_moving_average_lp_length];
    moving_average_lp_length_n7andback = moving_average_lp_length_n7andback(1:end-5);
    shuffled_MA_n7_Distribution = [shuffled_MA_n7_Distribution moving_average_lp_length_n7andback];
    %shuffled cumulative rewards
    shuffled_total_reward_indicator = total_reward_indicator;
    shuffled_total_reward_indicator = shuffled_total_reward_indicator(randperm(length(shuffled_total_reward_indicator)));
    shuffled_total_reward_indicator_dist = [shuffled_total_reward_indicator_dist shuffled_total_reward_indicator];    
     
    %independentyly shuffle ipi and he
    shuf_ipi1= ipi_nback_matrix(2:end,1); %remove the Nan at the start to avoid shuffling it
    shuffled_ipi1 = shuf_ipi1(randperm(length(shuf_ipi1)));
    %removed the original nan, now add it back in
    shuffled_ipi1 = [NaN; shuffled_ipi1];
    shuffled_ipi1_Distribution=[shuffled_ipi1_Distribution  shuffled_ipi1];  
    
      if isempty(ipi_nback_matrix(3:end,2)) %check to see if there are more than 1 lps
        shuffled_ipi2 = NaN;
        shuffled_ipi2_Distribution=[shuffled_ipi2_Distribution  shuffled_ipi2]; 
        shuffled_ipi3 = NaN;
        shuffled_ipi3_Distribution=[shuffled_ipi3_Distribution  shuffled_ipi3];
     else
            
    shuf_ipi2= ipi_nback_matrix(3:end,2);
    shuffled_ipi2 = shuf_ipi2(randperm(length(shuf_ipi2)));
    shuffled_ipi2 = [NaN; NaN; shuffled_ipi2];
    shuffled_ipi2_Distribution=[shuffled_ipi2_Distribution  shuffled_ipi2]; 
    %ipi n-3
    shuf_ipi3= ipi_nback_matrix(4:end,3);
    shuffled_ipi3 = shuf_ipi3(randperm(length(shuf_ipi3)));
    shuffled_ipi3 = [NaN; NaN; NaN; shuffled_ipi3];
    shuffled_ipi3_Distribution=[shuffled_ipi3_Distribution  shuffled_ipi3];  
    
    end
    %he n - 1          
    shuf_he1 = HE_n_1_Indicator(2:end);%again, remove Nan, shuffle, readd NaN
    shuf_he1 =  shuf_he1(randperm(length(shuf_he1)));
    shuf_he1 = [NaN; shuf_he1];
    shuffled_he1_Distribution =[shuffled_he1_Distribution shuf_he1];
    %he n -2
       %check for there only being a single LP
     if isempty( HE_n_2_Indicator(3:end))
    shuf_he2 = NaN;
    shuffled_he2_Distribution =[shuffled_he2_Distribution shuf_he2];
    shuf_he3 = NaN;
    shuffled_he3_Distribution =[shuffled_he3_Distribution shuf_he3];
    shuf_he4 = NaN;
    shuffled_he4_Distribution =[shuffled_he4_Distribution shuf_he4];

     else
    shuf_he2 = HE_n_2_Indicator(3:end);
    shuf_he2 =  shuf_he2(randperm(length(shuf_he2)));
    shuf_he2 = [NaN; NaN; shuf_he2];
    shuffled_he2_Distribution =[shuffled_he2_Distribution shuf_he2];
    %he n -3
    shuf_he3 = HE_n_3_Indicator(4:end);
    shuf_he3 =  shuf_he3(randperm(length(shuf_he3)));
    shuf_he3 = [NaN; NaN; NaN; shuf_he3];
    shuffled_he3_Distribution =[shuffled_he3_Distribution shuf_he3];
    %he n - 4
    shuf_he4 = HE_n_4_Indicator(5:end);
    shuf_he4 =  shuf_he4(randperm(length(shuf_he4)));
    shuf_he4 = [NaN; NaN; NaN; NaN; shuf_he4];
    shuffled_he4_Distribution =[shuffled_he4_Distribution shuf_he4];
     end
    %up-state
    Shuffled_Up_State_idx_n1 = Up_State_idx_n1(2:end);
    Shuffled_Up_State_idx_n1 = Shuffled_Up_State_idx_n1(randperm(length(Shuffled_Up_State_idx_n1)));
    Shuffled_Up_State_idx_n1 = [NaN; Shuffled_Up_State_idx_n1];
    Shuffled_Up_State_idx_n1_Distribution = [Shuffled_Up_State_idx_n1_Distribution Shuffled_Up_State_idx_n1];
    %stim indicators
    Shuffled_stim_ind_all = stim_ind_all(randperm(length(stim_ind_all)));
    if isempty(stim_ind_all_n1(2:end)) %checl for more than one lp
            Shuffled_stim_ind_all_n1 = [NaN];
            Shuffled_stim_ind_all_n2 = [NaN];
            Shuffled_stim_ind_all_n3 = [NaN];
            Shuffled_stim_ind_all_n4 = [NaN];
            Shuffled_stim_ind_all_n5 = [NaN];
            Shuffled_stim_ind_all_n6 = [NaN];

    else
    Shuffled_stim_ind_all_n1 = stim_ind_all_n1(2:end);
    Shuffled_stim_ind_all_n1 = Shuffled_stim_ind_all_n1(randperm(length(Shuffled_stim_ind_all_n1)));
    Shuffled_stim_ind_all_n1 = [NaN;Shuffled_stim_ind_all_n1];
    Shuffled_stim_ind_all_n2 = stim_ind_all_n2(3:end);
    Shuffled_stim_ind_all_n2 = Shuffled_stim_ind_all_n2(randperm(length(Shuffled_stim_ind_all_n2)));
    Shuffled_stim_ind_all_n2 = [NaN; NaN; Shuffled_stim_ind_all_n2];
    Shuffled_stim_ind_all_n3 = stim_ind_all_n3(4:end);
    Shuffled_stim_ind_all_n3 = Shuffled_stim_ind_all_n3(randperm(length(Shuffled_stim_ind_all_n3)));
    Shuffled_stim_ind_all_n3 = [NaN; NaN; NaN; Shuffled_stim_ind_all_n3];
    Shuffled_stim_ind_all_n4 = stim_ind_all_n4(5:end);
    Shuffled_stim_ind_all_n4 = Shuffled_stim_ind_all_n4(randperm(length(Shuffled_stim_ind_all_n4)));
    Shuffled_stim_ind_all_n4 = [NaN; NaN; NaN; NaN; Shuffled_stim_ind_all_n4];
    Shuffled_stim_ind_all_n5 = stim_ind_all_n5(6:end);
    Shuffled_stim_ind_all_n5 = Shuffled_stim_ind_all_n5(randperm(length(Shuffled_stim_ind_all_n5)));
    Shuffled_stim_ind_all_n5 = [NaN; NaN; NaN; NaN; NaN; Shuffled_stim_ind_all_n5];
            %check for if there are fewer than 5 values, due to low amounts of LPing then remove number of nans needed 
        if isempty(stim_ind_all_n5(6:end))   
    Shuffled_stim_ind_all_n5 = NaN(length(Lever_Press_Lengths),1);
    end
  
    Shuffled_stim_ind_all_n6 = stim_ind_all_n6(7:end);
    Shuffled_stim_ind_all_n6 = Shuffled_stim_ind_all_n6(randperm(length(Shuffled_stim_ind_all_n6)));
    Shuffled_stim_ind_all_n6 = [NaN; NaN; NaN; NaN; NaN; NaN; Shuffled_stim_ind_all_n6];
        %check for if there are fewer than 6 values, then remove number of nans needed 
    if isempty(stim_ind_all_n6(7:end))   
    Shuffled_stim_ind_all_n6 = NaN(length(Lever_Press_Lengths),1);
    end
    
    end %1 lp check end
    
    Shuffled_stim_ind_all_Distribution =[Shuffled_stim_ind_all_Distribution Shuffled_stim_ind_all]; 
    Shuffled_stim_ind_all_n1_Distribution =[Shuffled_stim_ind_all_n1_Distribution Shuffled_stim_ind_all_n1];    
    Shuffled_stim_ind_all_n2_Distribution =[Shuffled_stim_ind_all_n2_Distribution Shuffled_stim_ind_all_n2]; 
    Shuffled_stim_ind_all_n3_Distribution =[Shuffled_stim_ind_all_n3_Distribution Shuffled_stim_ind_all_n3];
    Shuffled_stim_ind_all_n4_Distribution =[Shuffled_stim_ind_all_n4_Distribution Shuffled_stim_ind_all_n4];
    Shuffled_stim_ind_all_n5_Distribution =[Shuffled_stim_ind_all_n5_Distribution Shuffled_stim_ind_all_n5];
    Shuffled_stim_ind_all_n6_Distribution =[Shuffled_stim_ind_all_n6_Distribution Shuffled_stim_ind_all_n6];

       %calulate shuffled durations after success and failure duration
       %changes. This will see how much change there is after
       %success/failure based purely on the distribution of durations,
       %rather than the order of presses
     Shuffled_Lever_Press_Rewarded_idx = find(Shuffled_Durs >= Data.Day(day).Mouse(mouse).Hold_Criteria); % Index of when successful Lever Press Ocurred
     Shuffled_Lever_Press_Rewarded_Lengths = Shuffled_Durs(Shuffled_Lever_Press_Rewarded_idx); % Length successful Lever Press
     Shuffled_Lever_Press_Not_Rewarded_idx = find(~(Shuffled_Durs >= Data.Day(day).Mouse(mouse).Hold_Criteria)); % Index of when unsuccessful Lever Press Ocurred
     Shuffled_Lever_Press_Not_Rewarded_Lengths = Shuffled_Durs(Shuffled_Lever_Press_Not_Rewarded_idx); % Length unsuccessful Lever Press
     Shuffled_Lever_Press_Length_After_Success_idx = Shuffled_Lever_Press_Rewarded_idx + 1; % Look for lever press index after success
     Shuffled_Lever_Press_Length_After_Failure_idx = Shuffled_Lever_Press_Not_Rewarded_idx + 1; % Look for lever press index after failure
      %for each of these, remove the last value if it is larger than the length of the original array
      %i.e., there is no n + 1 because its the final press in a shuffled
      %session
      if isempty(Shuffled_Lever_Press_Length_After_Success_idx)
          Shuffled_Lever_Press_Length_After_Success_idx =[];
          
      elseif Shuffled_Lever_Press_Length_After_Success_idx(end) > length(Shuffled_Durs)
     Shuffled_Lever_Press_Length_After_Success_idx(end) = [];
     Shuffled_Lever_Press_Rewarded_idx(end) =[];
      end
       if Shuffled_Lever_Press_Length_After_Failure_idx(end) > length(Shuffled_Durs)
     Shuffled_Lever_Press_Length_After_Failure_idx(end) = [];
     Shuffled_Lever_Press_Not_Rewarded_idx(end)  = [];
      end
              
          
        Shuffled_Lever_Press_Length_After_Success = Shuffled_Durs(Shuffled_Lever_Press_Length_After_Success_idx); % Lever Press Length at n+1 after success
        Shuffled_Lever_Press_Length_After_Failure = Shuffled_Durs(Shuffled_Lever_Press_Length_After_Failure_idx); % Lever Press Length at n+1 after failure
        Shuffled_Lever_Press_Length_After_Success_Change = Shuffled_Lever_Press_Length_After_Success - Shuffled_Durs(Shuffled_Lever_Press_Rewarded_idx);
        Shuffled_Lever_Press_Length_After_Failure_Change = Shuffled_Lever_Press_Length_After_Failure - Shuffled_Durs(Shuffled_Lever_Press_Not_Rewarded_idx);
%         
%calculate the mean duration after success/failure for each shuffle
mean_Shuffled_Lever_Press_Length_After_Success = mean(Shuffled_Lever_Press_Length_After_Success);
mean_Shuffled_Lever_Press_Length_After_Failure = mean(Shuffled_Lever_Press_Length_After_Failure);
std_Shuffled_Lever_Press_Length_After_Success_Dist = [std_Shuffled_Lever_Press_Length_After_Success_Dist std(Shuffled_Lever_Press_Length_After_Success)];
std_Shuffled_Lever_Press_Length_After_Failure_Dist = [std_Shuffled_Lever_Press_Length_After_Failure_Dist std(Shuffled_Lever_Press_Length_After_Failure)];   
mean_Shuffled_Lever_Press_Length_After_Success_Dist = [mean_Shuffled_Lever_Press_Length_After_Success_Dist mean_Shuffled_Lever_Press_Length_After_Success];
mean_Shuffled_Lever_Press_Length_After_Failure_Dist = [mean_Shuffled_Lever_Press_Length_After_Failure_Dist mean_Shuffled_Lever_Press_Length_After_Failure];         
mean_Shuffled_Lever_Press_Change_After_Success_Dist = [mean_Shuffled_Lever_Press_Change_After_Success_Dist mean(Shuffled_Lever_Press_Length_After_Success_Change)];
mean_Shuffled_Lever_Press_Change_After_Failure_Dist = [mean_Shuffled_Lever_Press_Change_After_Failure_Dist mean(Shuffled_Lever_Press_Length_After_Failure_Change)];         

 %use the shuffled lps to calculate up-state stuff. Then get distributions of proportion above and mean conseq.
 
 if isempty(Shuffled_Durs) %check to make sure there are more than 0 lps
    Up_State_idx_shuf = NaN;
    ilower_shuf = NaN;
    uppersum_shuf = NaN;
    lowersum_shuf = NaN;
 else
     
    if std(Shuffled_Durs) == 0 %if std is 0 because of only 1 LP, just code it as 1ms to avoid error
     [Up_State_idx_shuf,ilower_shuf,uppersum_shuf,lowersum_shuf] = cusum(Shuffled_Durs,2,1,mean(Shuffled_Durs),1,'all');
     
    else
     
     [Up_State_idx_shuf,ilower_shuf,uppersum_shuf,lowersum_shuf] = cusum(Shuffled_Durs,2,1,mean(Shuffled_Durs),std(Shuffled_Durs),'all');
    end  
 end
          %get the proportion of up state to the all LPS in the shuffled
          %version
         up_state_prop_shuf = length(Up_State_idx_shuf)/length(Shuffled_Durs);
       Up_State_Prop_Above_Shuf = [Up_State_Prop_Above_Shuf up_state_prop_shuf];
  %loop through the upstate  index to find the avg consecutive number of up state presses
sequence_above_shuf = 0;
conseq_shuf =0;
sequence_counter_shuf =[0];
for upper_it_shuf = 1:length(Up_State_idx_shuf)-1
   current_shuf = Up_State_idx_shuf(upper_it_shuf);
   if Up_State_idx_shuf(upper_it_shuf+1) == current_shuf +1
    
     conseq_shuf = conseq_shuf + 1;
     if conseq_shuf == 1
        sequence_above_shuf = sequence_above_shuf + 2;
       
     else 
        sequence_above_shuf = sequence_above_shuf + 1;

         continue
     end
     
   else
          sequence_counter_shuf = [sequence_counter_shuf conseq_shuf+1];
          conseq_shuf = 0;
        
       continue
   end        

end
%add last idx on since cant loop past the length
%remove the first 0
 sequence_counter_shuf = [sequence_counter_shuf conseq_shuf+1];
  sequence_counter_shuf(sequence_counter_shuf==0) =[];
  
sequence_counter_shuf(1) = sequence_counter_shuf(1) +1;
sequence_length_mean_including_ones_shuf = mean(sequence_counter_shuf);
Up_State_Consecutive_Shuf_Dist = [Up_State_Consecutive_Shuf_Dist sequence_counter_shuf];
Up_State_Consecutive_Shuf_Mean = [Up_State_Consecutive_Shuf_Mean sequence_length_mean_including_ones_shuf];
        
end
Data.Day(day).Mouse(mouse).shuffled_total_reward_indicator_dist = shuffled_total_reward_indicator_dist;
Data.Day(day).Mouse(mouse).shuffled_MA_n7_Distribution =shuffled_MA_n7_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_Distribution = Shuffled_stim_ind_all_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n1_Distribution = Shuffled_stim_ind_all_n1_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n2_Distribution = Shuffled_stim_ind_all_n2_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n3_Distribution = Shuffled_stim_ind_all_n3_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n4_Distribution = Shuffled_stim_ind_all_n4_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n5_Distribution = Shuffled_stim_ind_all_n5_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_stim_ind_all_n6_Distribution = Shuffled_stim_ind_all_n6_Distribution;

Data.Day(day).Mouse(mouse).shuffled_ipi1_Distribution = shuffled_ipi1_Distribution;
Data.Day(day).Mouse(mouse).shuffled_ipi2_Distribution = shuffled_ipi2_Distribution;
Data.Day(day).Mouse(mouse).shuffled_ipi3_Distribution = shuffled_ipi3_Distribution;
Data.Day(day).Mouse(mouse).shuffled_he1_Distribution = shuffled_he1_Distribution;
Data.Day(day).Mouse(mouse).shuffled_he2_Distribution = shuffled_he2_Distribution;
Data.Day(day).Mouse(mouse).shuffled_he3_Distribution = shuffled_he3_Distribution;
Data.Day(day).Mouse(mouse).shuffled_he4_Distribution = shuffled_he4_Distribution;
Data.Day(day).Mouse(mouse).shuffled_MA_Distribution = shuffled_MA_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_ts_Distribution = Shuffled_ts_Distribution;
Data.Day(day).Mouse(mouse).Shuffled_Up_State_idx_n1_Distribution = Shuffled_Up_State_idx_n1_Distribution;

Data.Day(day).Mouse(mouse).mean_Shuffled_Lever_Press_Length_After_Success_Dist = mean_Shuffled_Lever_Press_Length_After_Success_Dist;
Data.Day(day).Mouse(mouse).mean_Shuffled_Lever_Press_Length_After_Failure_Dist = mean_Shuffled_Lever_Press_Length_After_Failure_Dist;
Data.Day(day).Mouse(mouse).std_Shuffled_Lever_Press_Length_After_Success_Dist = std_Shuffled_Lever_Press_Length_After_Success_Dist;
Data.Day(day).Mouse(mouse).std_Shuffled_Lever_Press_Length_After_Failure_Dist = std_Shuffled_Lever_Press_Length_After_Failure_Dist;
Data.Day(day).Mouse(mouse).mean_Shuffled_Lever_Press_Change_After_Failure_Dist = mean_Shuffled_Lever_Press_Change_After_Failure_Dist;
Data.Day(day).Mouse(mouse).mean_Shuffled_Lever_Press_Change_After_Success_Dist = mean_Shuffled_Lever_Press_Change_After_Success_Dist;

Data.Day(day).Mouse(mouse).entropy_Shuffled_Lever_Press_Length_After_Success_Dist = entropy_Shuffled_Lever_Press_Length_After_Success_Dist;
Data.Day(day).Mouse(mouse).entropy_Shuffled_Lever_Press_Length_After_Failure_Dist = entropy_Shuffled_Lever_Press_Length_After_Failure_Dist;

Data.Day(day).Mouse(mouse).Up_State_Consecutive_Shuf_Dist = Up_State_Consecutive_Shuf_Dist;
Data.Day(day).Mouse(mouse).Up_State_Consecutive_Shuf_Mean = Up_State_Consecutive_Shuf_Mean;
Data.Day(day).Mouse(mouse).Up_State_Prop_Above_Shuf = Up_State_Prop_Above_Shuf;  

Up_State_Prop_Above_Shuf_Dist_Mean = mean(Up_State_Prop_Above_Shuf);
Up_State_Consecutive_Shuf_Mean_Dist_Mean = mean(Up_State_Consecutive_Shuf_Mean);
Data.Day(day).Mouse(mouse).Up_State_Prop_Above_Shuf_Dist_Mean =Up_State_Prop_Above_Shuf_Dist_Mean;
Data.Day(day).Mouse(mouse).Up_State_Consecutive_Shuf_Mean_Dist_Mean = Up_State_Consecutive_Shuf_Mean_Dist_Mean;

%now that we have created the shuffled lp orders, we need to generate
%n-backs for each shuffled LP order
shuf_n_back_Lengths_distribution = {};
shuf_Logical_n_back_Lengths_distribution = {};  
for shuffle_length = 1:size(Shuffled_Durs_Distribtuion,2)
    shuf_n_back_Lengths = [];
    shuf_n_lengths_array = [];
    shuf_Logical_n_lengths_array= [];
    shuf_Logical_n_back_Lengths =[];   
       
    current_shuffle_lps = Shuffled_Durs_Distribtuion(:,shuffle_length);
          shuf_press_indices = 1:length(current_shuffle_lps);
           shuf_press_indices = shuf_press_indices';
           
                
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
          shuf_Logical_n_lengths_array = shuf_n_lengths_array;          
     
         shuf_Logical_n_lengths_array(x:end) = shuf_Logical_n_lengths_array(x:end) >= Data.Day(day).Mouse(mouse).Hold_Criteria;
           
       end
       
% shuf_n_lengths_array
       shuf_n_back_Lengths = [shuf_n_back_Lengths; shuf_n_lengths_array]; 
      shuf_n_lengths_array =[];
    
      %the logical test converts all NaNs to 0, need to put them back so we
      %don't inappropriately put lever press prior to the start of the
      %session
               shuf_Logical_n_lengths_array = double(shuf_Logical_n_lengths_array);
               %if they make fewer than 10 lps, can sometimes be an issue
               
         shuf_Logical_n_back_Lengths = [shuf_Logical_n_back_Lengths; shuf_Logical_n_lengths_array];  
         shuf_Logical_n_lengths_array =[];
   
       end
       shuf_n_back_Lengths_distribution = [shuf_n_back_Lengths_distribution; shuf_n_back_Lengths];
       shuf_Logical_n_back_Lengths_distribution = [shuf_Logical_n_back_Lengths_distribution; shuf_Logical_n_back_Lengths];
    
end
    
Data.Day(day).Mouse(mouse).shuf_n_back_Lengths_distribution = shuf_n_back_Lengths_distribution;
Data.Day(day).Mouse(mouse).shuf_Logical_n_back_Lengths_distribution = shuf_Logical_n_back_Lengths_distribution;
Data.Day(day).Mouse(mouse).Shuffled_Durs_Distribtuion = Shuffled_Durs_Distribtuion;
      %% Convert LP lengths to logicals for successful or not
      Logical_LPs = Lever_Press_Lengths >= Data.Day(day).Mouse(mouse).Hold_Criteria;
      Data.Day(day).Mouse(mouse).Logical_LPs = Logical_LPs;

%% get cumulative rewards across a session
%the purpose of this section is to get performance across a session,
%grouped by cumulative rewards earned 
 %put the current lps (from one session, one animal) into matrix
%find the indicies of when reward occurred
Reward_indices = find(Logical_LPs);
%create variables that will get clear per animal/session    
cumulative_reward_lps_mouse =  {};
avg_duration_mouse =[ NaN(1,60)];
total_presses_mouse = [ NaN(1,60)];
correct_percent_mouse = [ NaN(1,60)];
std_duration_mouse = [NaN(1,60)];     
med_duration_mouse = [NaN(1,60)];       

%now use those indices to partition the lp matrix by rewards
%starting from index one
lp_index = 0;
for rew_it = 1:length(Reward_indices)
      %at each cumulative reward, calculate things like mean,med, std
      %duration of lever presses.
       lps_before_reward = Lever_Press_Lengths(lp_index+1:Reward_indices(rew_it));
       cumulative_reward_lps_mouse = [cumulative_reward_lps_mouse; lps_before_reward];
       avg_duration = mean(lps_before_reward);
       med_duration = median(lps_before_reward);
       total_presses = length(lps_before_reward);
       correct_percent = sum((lps_before_reward >= criteria_indicator(1)))/total_presses;
       std_presses = std(lps_before_reward);
       avg_duration_mouse(rew_it) = avg_duration;
       med_duration_mouse(rew_it) = med_duration;
       total_presses_mouse(rew_it) = total_presses;
       correct_percent_mouse(rew_it) =  correct_percent;
       std_duration_mouse(rew_it) = std_presses;
         %refresh the index to begin from +1 from the last reward, check to
         %make sure its not the last reward
         if Reward_indices(rew_it) == length(Lever_Press_Lengths)
             continue
         else
       lp_index = Reward_indices(rew_it);
         end
    end
 %after getting all the data for one animal, input the various matrices into cells  
  Data.Day(day).Mouse(mouse).avg_duration_by_reward = avg_duration_mouse;
  Data.Day(day).Mouse(mouse).total_presses_by_reward = total_presses_mouse;
  Data.Day(day).Mouse(mouse).correct_percent_by_reward = correct_percent_mouse;
  Data.Day(day).Mouse(mouse).cumulative_reward_lps_all_by_reward = cumulative_reward_lps_mouse;
  Data.Day(day).Mouse(mouse).std_duration_by_reward = std_duration_mouse;
  Data.Day(day).Mouse(mouse).med_duration_by_reward = med_duration_mouse;
     
%% sample shuffles
%code to shuffle the order of LPs,this one off shuffle is used for
%prediction comparisons 
Shuffled_Lever_Press_lengths = Lever_Press_Lengths;
Shuffled_Lever_Press_lengths =Shuffled_Lever_Press_lengths(randperm(length(Shuffled_Lever_Press_lengths)));
Data.Day(day).Mouse(mouse).Shuffled_Lever_Press_lengths = Shuffled_Lever_Press_lengths;

Shuffled_Logical_Lever_press = Logical_LPs;
Shuffled_Logical_Lever_press = Shuffled_Logical_Lever_press(randperm(length(Shuffled_Logical_Lever_press)));
Data.Day(day).Mouse(mouse).Shuffled_Logical_Lever_press = Shuffled_Logical_Lever_press;
       
%% get individual lp duration histograms per mouse
%100ms bins from 0 to 5000ms
[histogram_counts,histogram_edges,histogram_bin_indexs] = histcounts(Lever_Press_Lengths,'BinWidth',100,'BinLimits',[0 5000]); 
Data.Day(day).Mouse(mouse).histogram_counts = histogram_counts;
Data.Day(day).Mouse(mouse).histogram_bin_indexs = histogram_bin_indexs;
Data.Day(day).Mouse(mouse).histogram_edges = histogram_edges;

%% make an individual lme for each mouse/day
n_back_Lengths_mat = cell2mat(n_back_Lengths);
n_back_Lengths_mat = n_back_Lengths_mat';
%create a table that will include LP durations, n -1 through n - 10, and
%control covariates of cumulative rewards, time in session, and overall
%crieria%
T10_Logical_and_Continuous = [Lever_Press_Lengths n_back_Lengths_mat total_reward_indicator Lever_Press_ts crit_percent_indicator ];
if isempty(T10_Logical_and_Continuous);
emptycheck = 1;
elseif size(T10_Logical_and_Continuous,1) <10
    shortcheck = 1;
else
    
T10_Logical_and_Continuous=array2table(T10_Logical_and_Continuous,'VariableNames',{'LP_Durations_All','n_minus_one_Durations_All',...
    'n_minus_two_Durations_All', 'n_minus_three_Durations_All', 'n_minus_four_Durations_All',...
    'n_minus_five_Durations_All', 'n_minus_six_Durations_All', 'n_minus_seven_Durations_All',...
    'n_minus_eight_Durations_All', 'n_minus_nine_Durations_All', 'n_minus_ten_Durations_All',...
    'total_reward_indicator_All','LP_Timestamps_All','criteria_percent_indicator'});
%predict duration from n-back duration and session covariates. This is NOT
%mixed (ie no random effects), since we are only looking at a single
%session for a single mouse

model_spec_10_ma_totrew = 'LP_Durations_All~LP_Timestamps_All+n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All';
% model_spec_10_ma_totrew = 'LP_Durations_All~n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All ';
lme_10_ma_totrew = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew);
lme_10_ma_totrew_se = lme_10_ma_totrew.Coefficients.SE;
lme_10_ma_totrew_coef =  lme_10_ma_totrew.Coefficients.Estimate;
lme_10_ma_totrew_name =  lme_10_ma_totrew.Coefficients.Name;
lme_10_ma_totrew_pval =  lme_10_ma_totrew.Coefficients.pValue;

%predict durations based on the model. We want to see how good the model is
%at predicting behavior for individual mice
[yhat yhatCI yhatDF]= predict(lme_10_ma_totrew,'Simultaneous',false);
 correctish_pred =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI(:,1);
 correctCI_prop = sum(correctish_pred)/length(T10_Logical_and_Continuous.LP_Durations_All)  ;
%put this with %criteria in an array, just to see if performance correlates with predictability from the model  
correctCIprop_criteria = [correctCI_prop crit_percent_indicator(1)];

r_squared_criteria = [lme_10_ma_totrew.Rsquared.Ordinary crit_percent_indicator(1)/100];
r_squared_adjusted_criteria = [lme_10_ma_totrew.Rsquared.Adjusted crit_percent_indicator(1)/100];

%lets calculate the Coef. partial determination for n - 1 duration. We will
%eventually correlate this with criteria% to see if how much variance is
%explained by n -1 correlates with how well animals perform
model_spec_10_ma_totrew_minus_n1 = 'LP_Durations_All~LP_Timestamps_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All+n_minus_five_Durations_All+n_minus_six_Durations_All+n_minus_seven_Durations_All+n_minus_eight_Durations_All+n_minus_nine_Durations_All+n_minus_ten_Durations_All';
% model_spec_10_ma_totrew_minus_n1 = 'LP_Durations_All~ n_minus_two_Durations_All + n_minus_three_Durations_All ';
lme_10_ma_totrew_minus_n1 = fitlme(T10_Logical_and_Continuous,model_spec_10_ma_totrew_minus_n1);

indiv_CPD_n1_simple = lme_10_ma_totrew.Rsquared.Ordinary - lme_10_ma_totrew_minus_n1.Rsquared.Ordinary;
indiv_CPD_n1_criteria_simple = [indiv_CPD_n1_simple crit_percent_indicator(1)/100];

%add in some of the other possible events to the table
T10_Logical_and_Continuous.HE_n_2_Indicator_All = HE_n_2_Indicator;
T10_Logical_and_Continuous.HE_n_1_Indicator_All = HE_n_1_Indicator;
T10_Logical_and_Continuous.ipi1 = n_to_n_back_ipi_all{1,1}';
T10_Logical_and_Continuous.ipi2 = n_to_n_back_ipi_all{2,1}';;
T10_Logical_and_Continuous.n_minus_one_All = Logical_n_back_Lengths{1,1}';
T10_Logical_and_Continuous.moving_average_lp_length_n7andback_All = moving_average_lp_length_n7andback;
T10_Logical_and_Continuous.up_state_idx_n1_All = Up_State_idx_n1;
T10_Logical_and_Continuous.HE_n_1_Indicator_All = categorical(T10_Logical_and_Continuous.HE_n_1_Indicator_All );
T10_Logical_and_Continuous.n_minus_one_All = categorical(T10_Logical_and_Continuous.n_minus_one_All);

%build a model that includes these other terms from BIC selection of all
%data. Remove criteria% since that is the same for all lps on a given day
% model_spec_moderate_remove = 'LP_Durations_All~ n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All ';
model_spec_moderate_remove  = 'LP_Durations_All ~ moving_average_lp_length_n7andback_All:LP_Timestamps_All + moving_average_lp_length_n7andback_All:n_minus_one_All + moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All + ipi1:moving_average_lp_length_n7andback_All  + n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All +    LP_Timestamps_All ';
lme__moderate_remove = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove);
lme__moderate_remove_coef =  lme__moderate_remove.Coefficients.Estimate;
lme__moderate_remove_se = lme__moderate_remove.Coefficients.SE;

%predict durations given the more complex model
[yhat yhatCI yhatDF]= predict(lme__moderate_remove,'Simultaneous',false);
correctish_pred_complex =T10_Logical_and_Continuous.LP_Durations_All <=yhatCI(:,2) &  T10_Logical_and_Continuous.LP_Durations_All >=yhatCI(:,1);
correctCI_prop_complex = sum(correctish_pred_complex)/length(T10_Logical_and_Continuous.LP_Durations_All)  ;
correctCIprop_criteria_complex = [correctCI_prop_complex crit_percent_indicator(1)];
r_squared_complex_criteria = [lme__moderate_remove.Rsquared.Ordinary crit_percent_indicator(1)/100];
r_squared_adjusted_complex_criteria = [lme__moderate_remove.Rsquared.Adjusted crit_percent_indicator(1)/100];
 
%remove all N-1 duration terms to see how much n1 contributes overall
model_spec_moderate_remove_minus_n1 = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All +  moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All + moving_average_lp_length_n7andback_All:n_minus_one_All +   HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All +  n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All ';
% lme__moderate_remove_minus_n1 = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_minus_n1);
lme__moderate_remove_minus_n1 = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove);

indiv_CPD_n1_complex = lme__moderate_remove.Rsquared.Ordinary - lme__moderate_remove_minus_n1.Rsquared.Ordinary;
indiv_CPD_n1_criteria_complex = [indiv_CPD_n1_complex crit_percent_indicator(1)/100];
%ditto but for the moving average term
model_spec_moderate_remove_minus_MA = 'LP_Durations_All ~ n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All + n_minus_one_Durations_All:n_minus_one_All  +  HE_n_1_Indicator_All + ipi1 + n_minus_one_All  +  ipi2:n_minus_two_Durations_All +   ipi2  +  n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All ';
lme__moderate_remove_minus_MA = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_minus_MA);
% lme__moderate_remove_minus_MA = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove);

indiv_CPD_MA_complex = lme__moderate_remove.Rsquared.Ordinary - lme__moderate_remove_minus_MA.Rsquared.Ordinary;
indiv_CPD_MA_criteria_complex = [indiv_CPD_MA_complex crit_percent_indicator(1)/100];
%remove n-1 reward and interaction with n1 and MA
model_spec_moderate_remove_minus_n1rew = 'LP_Durations_All ~  moving_average_lp_length_n7andback_All:LP_Timestamps_All +  moving_average_lp_length_n7andback_All:HE_n_1_Indicator_All +  ipi1:moving_average_lp_length_n7andback_All +  n_minus_one_Durations_All:LP_Timestamps_All + n_minus_one_Durations_All:HE_n_1_Indicator_All +  ipi1:n_minus_one_Durations_All +  HE_n_1_Indicator_All + ipi1 +  ipi2:n_minus_two_Durations_All +   ipi2  + moving_average_lp_length_n7andback_All + n_minus_one_Durations_All + n_minus_two_Durations_All + n_minus_three_Durations_All + n_minus_four_Durations_All + n_minus_five_Durations_All + n_minus_six_Durations_All +   LP_Timestamps_All ';
lme__moderate_remove_minus_n1rew = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove_minus_n1rew);
% lme__moderate_remove_minus_n1rew = fitlme(T10_Logical_and_Continuous,model_spec_moderate_remove);

indiv_CPD_n1rew_complex = lme__moderate_remove.Rsquared.Ordinary - lme__moderate_remove_minus_n1rew.Rsquared.Ordinary;
indiv_CPD_n1rew_criteria_complex = [indiv_CPD_n1rew_complex crit_percent_indicator(1)/100];

%save predictions and the beta coef/se for the simple model
Data.Day(day).Mouse(mouse).correctCI_prop =correctCI_prop;
Data.Day(day).Mouse(mouse).correctCIprop_criteria =correctCIprop_criteria;
Data.Day(day).Mouse(mouse).lme_10_indiv_se =lme_10_ma_totrew_se;
Data.Day(day).Mouse(mouse).lme_10_indiv_coef =lme_10_ma_totrew_coef;
%save predictions and coefs of complex model
Data.Day(day).Mouse(mouse).lme__moderate_remove_coef =  lme__moderate_remove_coef ;
Data.Day(day).Mouse(mouse).lme__moderate_remove_se = lme__moderate_remove_se;

Data.Day(day).Mouse(mouse).correctCI_prop_complex =correctCI_prop_complex;
Data.Day(day).Mouse(mouse).correctCIprop_criteria_complex =correctCIprop_criteria_complex;
%save r^2 values reltaive to criteria prop: does model fit correlate with
%%correct?
Data.Day(day).Mouse(mouse).r_squared_criteria = r_squared_criteria;
Data.Day(day).Mouse(mouse).r_squared_adjusted_criteria = r_squared_adjusted_criteria;
Data.Day(day).Mouse(mouse).r_squared_complex_criteria = r_squared_complex_criteria;
Data.Day(day).Mouse(mouse).r_squared_adjusted_complex_criteria = r_squared_adjusted_complex_criteria;
%save cpd for n -1 dur in simple, complex, and mov avg in complex
Data.Day(day).Mouse(mouse).indiv_CPD_n1_criteria_simple =  indiv_CPD_n1_criteria_simple ;
Data.Day(day).Mouse(mouse).indiv_CPD_n1_criteria_complex =  indiv_CPD_n1_criteria_complex ;
Data.Day(day).Mouse(mouse).indiv_CPD_MA_criteria_complex =  indiv_CPD_MA_criteria_complex ;
Data.Day(day).Mouse(mouse).indiv_CPD_n1rew_criteria_complex =  indiv_CPD_n1rew_criteria_complex ;

end

%lets also save the total time in a session and the number of rewards
%earned to be used for correlations
total_time = time_with_events(end,2);
total_time_and_rewards = [total_time subject_vars(REIN_idx) ];
Data.Day(day).Mouse(mouse).total_time_and_rewards = total_time_and_rewards;

%% Find the Difference between adjacent LPS
LP_Diffs = diff(Lever_Press_Lengths);
LP_Diffs = LP_Diffs./Lever_Press_Lengths(1:end-1);

med_LP_diff = median(LP_Diffs);
mean_LP_diff = mean(LP_Diffs);
std_LP_diff = std(LP_Diffs);

Data.Day(day).Mouse(mouse).LP_diffs = LP_Diffs;
Data.Day(day).Mouse(mouse).mean_LP_diff = mean_LP_diff;
Data.Day(day).Mouse(mouse).std_LP_diff = std_LP_diff;
Data.Day(day).Mouse(mouse).med_LP_diff = med_LP_diff;

LP_Increase = find(LP_Diffs >0);
LP_Decrease = find(LP_Diffs <0);

inc_prop = length(LP_Increase)/length(LP_Diffs);
dec_prop = length(LP_Decrease)/length(LP_Diffs);
dec_diff_mean = median(LP_Diffs(LP_Decrease));
inc_diff_mean = median(LP_Diffs(LP_Increase));
Data.Day(day).Mouse(mouse).inc_prop = inc_prop;
Data.Day(day).Mouse(mouse).dec_prop = dec_prop;
Data.Day(day).Mouse(mouse).dec_diff_mean = dec_diff_mean;
Data.Day(day).Mouse(mouse).inc_diff_mean = inc_diff_mean;
  
    end %end of the file loop
    
end %end of the folder loop
    cd(PathName_Folders)

%% ----- Choose where you save your extracted MEDPC data -----%%
cd(PathName_Master) % Change directory to one above each day's folder
save(Experiment_name , 'Data', '-v7.3') %need to use old 7.3 format for large file sizes. Becomes an issue with large sample sizes/shuffles
toc%timer
end
