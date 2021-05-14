%% This is the master script for GCaMP photometry analysis in Hold Down
%Last updated 5/10/21 by Drew Schreiner

%This script will run other functions to extract and analyze the data
%First we set some variables below, such as the criteria in ms, and the day
%of training. Then we point to the folder in which we have the csv files
%for photometry
%These files should include: a behavior csv (all behavioral data and
%timestamps) and a photometry csv (all photometry data and timestamps for
%each mouse/session
%This master script will loop through each individual's pair of
%behavior/photometry data
%You need to make sure that the behavior file ends in "analog" while the
%photometry end in "photometry", as matlab will search for csvs with those
%endings to extract data from.

%Functions Run by master_script:
%[GCAMP] = Photometry_r(GCAMP.mouseID, GCAMP.training_day);
% [photometry_data] = photometry_extract_r(CSVFiles.name);

    %Extracts photometry signal
% GCAMP.beh_data = behavior_extract_r(CSVFiles.name);
    %Extracts Behavioral Data
% [GCAMP] = GCAMP_data_extract_r(GCAMP);
    %This aligns the extraced photometry and behavioral data to get
    %behaviorally relevant photometry information (e.g., activity during
    %lever pressing
% [GCAMP] = GCAMP_plot_with_baseline_r(GCAMP, base_time_start, base_time_end, time_end);
    %These three will plot individual subject/day data, calculating
    %baselines, aligning activity to lever pressing, and make plots based
    %on whether presses were successful (exceeded criteria) or not.    
% [GCAMP] = GCAMP_Linear_Regressions_r(GCAMP, base_time_start, base_time_end, time_end);
    %This function will generate individual subject data in a table to be
    %used for linear mixed effects models relating behavior to gcamp
    %activity it specific timepoints.
 
 %For these two, you only need the grouped_gcamp data structure   
 % GCAMP_plot_with_baseline_Z_scores_r(Grouped_GCAMP, base_time_start, base_time_end,time_end)
    %This function will take as inpur the Grouped_GCAMP structure, which
    %will have as many mice/days together as were run in an individual
    %folder. From here, it will put all this data together and graph it   
% GCAMP_grand_regression_indivshuffles_r(Grouped_GCAMP, base_time_start, base_time_end,time_end)
    %This function takes the individual subject linear regression data
    %(tables) and puts it all together into one large table with indicator
    %variables for mouse/session. This is used to generate Linear Mixed
    %Effect models to relate activity at various timepoints to current and
    %prior behavior.
    
%% Organization info
tic
clear all
Criteria = 800; % in ms
%name the day of training, direct to dir of files/save dir
training_day = 'cie-800ms1';
PathName_Folder = 'C:\Users\dcschrei\Desktop\62-1 M2-DMS CIE GCAMP\0800MS1\cie'
cd(PathName_Folder);
GCAMP_Save_Dir = 'C:\Users\dcschrei\Desktop\62-1 M2-DMS CIE GCAMP\0800MS1\cie'
indiv_files = dir('*analog.csv');

%loop through the number of analog files, to get number of animal/day
mice ={};
for i = 1:length(indiv_files)
   indiv_mouse = indiv_files(i).name(1:17);
   mice = [mice indiv_mouse];
end
%for each photometry/behavioral pair
for i = 1:length(mice)
%Create the GCAMP data structure. This will hold all the data by mouse
GCAMP.mouseID = mice{i};
GCAMP.training_day = training_day;
GCAMP.Criteria = Criteria;
%This takes the first 4 digits of the filename to assign a mouse id
GCAMP.onlyID = GCAMP.mouseID(1:4);
%% Extract photometry and behavior data from .csv files
[GCAMP] = Photometry_r(GCAMP.mouseID, GCAMP.training_day);
%Find the mouse's analog behavioral data
CSVFiles = dir([GCAMP.mouseID '*analog.csv']);
%Extract Behavioral data from the analog csv file 
%Gets timestamps, Lever press, Headentry, and Reward
GCAMP.beh_data = behavior_extract_r(CSVFiles.name);

%Save the Behavioral and photometry data in the GCAMP structure in the save directory
cd(GCAMP_Save_Dir)
save(['GCAMP_' GCAMP.mouseID '_' GCAMP.training_day], 'GCAMP');
%Go pack to the folder with the data files
cd(PathName_Folder);
%Repeat for each mouse
end
%% Save data structure (100% Reward Probability)
%Go back to GCAMP_Save_Dir to open those files and run analysis
cd(GCAMP_Save_Dir)

for i = 1:length(mice)
    load(['GCAMP_' mice{i} '_' GCAMP.training_day])    
GCAMP.training_day = training_day;
GCAMP.Criteria = Criteria;
GCAMP.onlyID = GCAMP.mouseID(1:4);
%This function analyzes the raw behavioral data
[GCAMP] = GCAMP_data_extract_r(GCAMP);
%% Plot perievent (lever press) GCAMP data (WITH BASELINE)
%optional line to suppress plotting
% set(groot,'defaultFigureVisible','off')

%Set the window for the baseline and perievent end
%currently baseline is -15s to -5s, and event is -5s to +5s
base_time_start = -15;
base_time_end = -5;
time_end = 5;
%Plot perievent data. Align GCAMP signal with beh_data. Transform GCAMP
%signal (e.g., use baseline activity to z-score)
[GCAMP] = GCAMP_plot_with_baseline_r(GCAMP, base_time_start, base_time_end, time_end);

%Take the gcamp and beh data and build linear mixed effect models
%mostly to predict activity given behavior, and then save into table formt
%to create grand LMEs with all subjects/days of interest
[GCAMP] = GCAMP_Linear_Regressions_r(GCAMP, base_time_start, base_time_end, time_end);

%% Save individual GCAMP file and add to grouped GCAMP data structure
cd(GCAMP_Save_Dir)

save(['GCAMP_' GCAMP.mouseID '_' GCAMP.training_day], 'GCAMP');
% close all
%Add individual mouse GCAMP structure to a Grouped_GCAMP structure that
%will include all mice
Grouped_GCAMP.Mice{i} = GCAMP;

end
%% Plot Grouped Z-Scores
%Takes Grouped_GCAMP and plops mouse data together for summary graphs
GCAMP_plot_with_baseline_Z_scores_r(Grouped_GCAMP, base_time_start, base_time_end,time_end)

%change to save directory and save grouped data structure, use old Matlab format for large file
%size
directory = GCAMP_Save_Dir;
% save the grouped_gcamp structure
save(['GCAMP_Grouped-' Grouped_GCAMP.Mice{1}.training_day], 'Grouped_GCAMP','-v7.3');

%run the regression models using ALL the mouse data
GCAMP_grand_regression_indivshuffles_r(Grouped_GCAMP)


toc
