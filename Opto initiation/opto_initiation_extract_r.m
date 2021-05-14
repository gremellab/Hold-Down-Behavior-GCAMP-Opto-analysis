function MS8RETRAINanalogOPTO = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  MS8RETRAINANALOGOPTO = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the numeric data.
%
%  MS8RETRAINANALOGOPTO = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  MS8RETRAINanalogOPTO = importfile("C:\Users\dcschrei\Desktop\video-tracking optoinitiation\67-1 3412 day 2 test\3412-1600MS8-RETRAIN-analog-OPTO.csv", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 05-Jun-2020 08:19:32

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Item1", "Item2"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
MS8RETRAINanalogOPTO = readtable(filename, opts);

%% Convert to output type
MS8RETRAINanalogOPTO = table2array(MS8RETRAINanalogOPTO);
end