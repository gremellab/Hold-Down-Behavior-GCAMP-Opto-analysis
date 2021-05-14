%This master script will run functions to analyze the opto initiation data
%ultimately we want to build LMEs to see if light affects the contribution
%of prior experience

%Takes as input 4 excel CSVs generated during the task:
%...*analog file, which has the behavioral data and timestamps
%...*tracking file, which has data from the overhead camera position
%tracking
%...*opto file, which has timestamps and binary for when stimulation occurred (since it is triggered 50% of the time, 0's indicate when conditions where met but no stim occurred due to chance)
%...*zone file, which has timestamps for when mice entered and left the
%lever press initiation zone

%Then, extracts the data from these csvs using matlab generated fxs:
%beh_extract_r
% vid_extract_r
% opto_initiation_extract_r
% opto_zone_extract_z

%From there, this data is brought together and analyzed using 
%beh_and_vid_data_analysis_r

%Then this information is saved, and all mice are analzyed together via
%Grouped_Vid_Data_Function_r
%last Updated 5/11/21 Drew Schreiner

%% Get the data imported 
tic
clear all

Criteria = 1600; % in ms
training_day = 'opto-1600m1-6-alllatencies-revisedcode-alllatency';

PathName_Folder = 'C:\Users\dcschrei\Desktop\video-tracking optoinitiation\67-1-67-2 together\all stim no crash correct days\Opto Revised Code';
Data_Save_Dir =   'C:\Users\dcschrei\Desktop\video-tracking optoinitiation\67-1-67-2 together\all stim no crash correct days\Opto Revised Code';
indiv_files = dir('*ANALOG.csv');

%loop through the number of analog files, to get number of animal/day
mice ={};
for i = 1:length(indiv_files)
   indiv_mouse = indiv_files(i).name(1:17);
   mice = [mice indiv_mouse];
end
                                                                                                                                                                                                                                                                                                                  
%Import the 4 csvs for each mouse and day
for i = 1:length(mice)
Vid_Data.mouseID = mice{i};
Vid_Data.training_day = training_day;
Vid_Data.day_indicator = mice{i}(12:13);
Vid_Data.Criteria = Criteria;
Vid_Data.onlyID = Vid_Data.mouseID(1:4);

CSVFiles_analog = dir([Vid_Data.mouseID '*analog.csv']);
beh_data_raw = beh_extract_r(CSVFiles_analog.name);
Vid_Data.beh_data_raw = beh_data_raw;

CSVFiles_track = dir([Vid_Data.mouseID '*tracking.csv']);
vid_data_raw = vid_extract_r(CSVFiles_track.name);
Vid_Data.vid_data_raw = vid_data_raw;

CSVFiles_opto = dir([Vid_Data.mouseID '*OPTO.csv']);
opto_data_raw = opto_initiation_extract_r(CSVFiles_opto.name);
Vid_Data.opto_data_raw = opto_data_raw;

CSVFiles_zone = dir([Vid_Data.mouseID '*ZONE.csv']);
zone_data_raw = opto_zone_extract_r(CSVFiles_zone.name);
Vid_Data.zone_data_raw = zone_data_raw;

%now the data is all imported, save into Vid_Data, data structure
save(['Vid_Data_' Vid_Data.mouseID '_' Vid_Data.training_day], 'Vid_Data');
cd(PathName_Folder);
end


%Now that the data is imported, reload each Vid_Data structure and run
%analysis on it
for i = 1:length(mice)
load(['Vid_Data_' mice{i} '_' Vid_Data.training_day])

Vid_Data.Criteria = Criteria;
Vid_Data.mouseID = mice{i};
Vid_Data.training_day = training_day;
Vid_Data.day_indicator = mice{i}(12:13);

% run the vid and beh analysis function
[Vid_Data] = beh_and_vid_data_analysis_r(Vid_Data.beh_data_raw, Vid_Data.vid_data_raw, Vid_Data.opto_data_raw, Vid_Data.zone_data_raw, Vid_Data.Criteria, Vid_Data.mouseID, Vid_Data.training_day);

% save the vid_data structure into a folder for your experiment
Vid_Data.training_day = training_day;
Vid_Data.Criteria = Criteria;
Vid_Data.mouseID = mice{i};
Vid_Data.day_indicator = mice{i}(12:13);

save(['Vid_Data_' Vid_Data.mouseID '_' Vid_Data.training_day], 'Vid_Data');
%save into a grouped data strucutre, which will have all mice
Grouped_Vid_Data.Mice{i} = Vid_Data;

end

save(['Grouped_Vid_Data-' Grouped_Vid_Data.Mice{1}.training_day], 'Grouped_Vid_Data','-v7.3');

%analyze the grouped data 
Grouped_Vid_Data_Function_r(Grouped_Vid_Data)

toc
