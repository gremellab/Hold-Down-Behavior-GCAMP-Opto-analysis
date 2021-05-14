function [GCAMP] = Photometry_r(mouseID,training_day)

%extracts the photometry data from an excel csv
%if using isobestic, finds which sample is 470 and 405/410
%and fits double exp for decay
CSVFiles = dir([mouseID '*photometry.csv']);

%automated matlab fx to extract columns from excel csv
[photometry_data] = photometry_extract_r(CSVFiles.name);

% name gcamp data wave and timestamps wavetime 
% delete first row if data has headers
raw_gcamptimestamps = photometry_data(:,1);
raw_gcampdata = photometry_data(:,2);
% if there are an odd # of rows after removing headers, delete the last row
% in both wave and wavetime
A = raw_gcampdata(1:2:end,:);
B = raw_gcampdata(2:2:end,:);

if length(A) > length(B)
  A = A(1:end-1,:);
elseif length(A) < length(B)
  B = B(1:end-1,:);
end

if sum(A>B) > (length(A)*.99)
  LED1 = A;
  LED2 = B;
else
  LED1 = B;
  LED2 = A;
end


% Timestamps for gcampdata
if LED1 == A
    k = 1;
elseif LED1 == B
    k = 2;
end
gcampdata_timestamps = raw_gcamptimestamps(k:2:end,:);

%% Double exp fit
SR = 20;
time = (1:length(LED1))/SR;
Y_exp_fit_all = fit(time',LED1,'exp2');
f_0 = Y_exp_fit_all(0);
normDat = (LED1 * f_0) ./ Y_exp_fit_all(time);

%if using 405 isobestic
  % Normalize to fit for 405
%   gcampdata = (LED1-fitLED)./fitLED;
% normDat = (dat1 - controlFit)./ controlFit; %this gives deltaF/F
% normDatPer = normDat * 100; % get %
gcampdata = normDat;

%% Save data to GCAMP structure
GCAMP.mouseID = mouseID;
GCAMP.training_day = training_day;
GCAMP.gcampdata = gcampdata;
GCAMP.gcampdata_timestamps = gcampdata_timestamps;
GCAMP.LED1 = LED1;
GCAMP.LED2 = LED2;
%GCAMP.fitLED = fitLED;
GCAMP.fitLED = Y_exp_fit_all;%controlFit;
GCAMP.raw_gcamptimestamps = raw_gcamptimestamps;
GCAMP.raw_gcampdata = raw_gcampdata;
GCAMP.SR = SR;
%% Plot raw GCAMP data to look for artifacts (e.g., sudden large drop indicative of de-coupling)

red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
cyan = [0.3010, 0.7450, 0.9330];
gray1 = [.7 .7 .7];
gray2 = [.8 .8 .8];

figure('Position',[100, 100, 800, 400])
hold on;
p1 = plot(time, GCAMP.LED1,'color',red,'LineWidth',2);
p2 = plot(time, gcampdata,'color',cyan,'LineWidth',2);
p3 = plot(time, Y_exp_fit_all(time),'color',green,'LineWidth',10);
title('Exponential Fit vs Raw','fontsize',16);
ylabel('F','fontsize',16);
axis tight;
legend([p1 p2 p3], {'Raw 470', 'New Normalized 470','Exp Fitted 470'});

end

