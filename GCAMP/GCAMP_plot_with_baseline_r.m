function [GCAMP] = GCAMP_plot_with_baseline_r(GCAMP, base_time_start, base_time_end, time_end)
%GCAMP_PLOT Plots perievent (lever press) GCAMP data relative TO BASELINE
%last updated 5/10/21 by Drew Schreiner
%10th% running window, used for checking df/f change
win = GCAMP.SR * 15;
p = 10;
[ y ] = running_percentile_filter(GCAMP.gcampdata,win,0,p);
delta_F = GCAMP.gcampdata - y;
GCAMP.delta_F = (delta_F ./ y) .* 100;
GCAMP.delta_F_z = zscore(delta_F);

%check to see that there is at least 1% change from baseline at 97.5
%prctile (as in Markowitz)
one_percent_check = prctile(GCAMP.delta_F ,97.5);
if one_percent_check < 1
   GCAMP.check = 'bad'
else
    GCAMP.check = 'good'   
end

%% Express data as baseline change
% you'll be putting each trial into a matrix with dimensions [trial * time]
raw_baseline_mean = [];
raw_baseline_std = [];
Closest_LP_ON_dist_all =[];
raw_F_LP_ON =[];
raw_F_LP_OFF =[];
off_raw_baseline_mean = [];
off_raw_baseline_std = [];
delta_F_LP_OFF_unsmoothed =[];
All_LP_ON_data_Baselines_raw=[];

%get event aligned data, go to a baseline period (-15 to -5s prior to
%event%(lp onset or offset)), take the mean and SD of baseline and use to z-score
%the rest of the trace from -5 to +5s relative to the event.
for i = 1:length(GCAMP.LP_ON_timestamps) %for each LP_ON timestamp
    % Find index of this LP_ON & OFF timestamp in the vector of all of your
    % timestamps
    [Closest_LP_ON_idx,Closest_LP_ON_dist] = nearestpoint(GCAMP.LP_ON_timestamps(i),GCAMP.gcampdata_timestamps);
    Closest_LP_OFF_idx = nearestpoint(GCAMP.LP_OFF_timestamps(i),GCAMP.gcampdata_timestamps);
    Closest_LP_ON_dist_all = [Closest_LP_ON_dist_all; Closest_LP_ON_dist];

    raw_F_LP_ON_individual_press = GCAMP.gcampdata(Closest_LP_ON_idx + base_time_end * GCAMP.SR : Closest_LP_ON_idx + time_end * GCAMP.SR)';
    raw_F_LP_OFF_individual_press = GCAMP.gcampdata(Closest_LP_OFF_idx + base_time_end * GCAMP.SR : Closest_LP_OFF_idx + time_end * GCAMP.SR)';
    raw_F_LP_ON = [raw_F_LP_ON; raw_F_LP_ON_individual_press];
    raw_F_LP_OFF = [raw_F_LP_OFF; raw_F_LP_OFF_individual_press];
               
    raw_baseline =  GCAMP.gcampdata(Closest_LP_ON_idx + base_time_start * GCAMP.SR : Closest_LP_ON_idx + base_time_end * GCAMP.SR);
    All_LP_ON_data_Baselines_raw = [All_LP_ON_data_Baselines_raw ; raw_baseline'];
    raw_baseline_mean = [raw_baseline_mean; mean(raw_baseline)];
    raw_baseline_std = [raw_baseline_std; std(raw_baseline)];
    raw_baseline_off =  GCAMP.gcampdata(Closest_LP_OFF_idx + base_time_start * GCAMP.SR : Closest_LP_OFF_idx + base_time_end * GCAMP.SR);
    off_raw_baseline_mean = [off_raw_baseline_mean; mean(raw_baseline_off)];
    off_raw_baseline_std = [off_raw_baseline_std; std(raw_baseline_off)];
end

%align to reinforcers
raw_F_RE_ON =[];
for i = 1:length(GCAMP.RE_ON_timestamps)
    Closest_RE_ON_idx = nearestpoint(GCAMP.RE_ON_timestamps(i),GCAMP.gcampdata_timestamps);
    raw_F_RE_ON_individual = GCAMP.gcampdata(Closest_RE_ON_idx + base_time_end * GCAMP.SR : Closest_RE_ON_idx + time_end * GCAMP.SR)';
    raw_F_RE_ON =[raw_F_RE_ON; raw_F_RE_ON_individual];
end
GCAMP.raw_F_RE_ON = raw_F_RE_ON;

%% z-score to baseline stuff
%the one without raw = df/f, raw is raw signal z-scored to baseline
baseline_z_score_LP_On_Raw = (raw_F_LP_ON -raw_baseline_mean)./raw_baseline_std;
raw_baseline_z_score_Met = baseline_z_score_LP_On_Raw(GCAMP.Criteria_met, :);
raw_baseline_z_score_Fail = baseline_z_score_LP_On_Raw(GCAMP.Criteria_fail, :);

%for the offset, we are actually going to use the baseline for lp onset.
%That way we will use the same baseline for onset/offset for a given lp
baseline_z_score_LP_OFF_Raw = (raw_F_LP_OFF -raw_baseline_mean)./raw_baseline_std;
baseline_z_score_LP_OFF_Raw_Met = baseline_z_score_LP_OFF_Raw(GCAMP.Criteria_met,:);
baseline_z_score_LP_OFF_Raw_Fail = baseline_z_score_LP_OFF_Raw(GCAMP.Criteria_fail,:);

%% Reinforced Lever Presses (Criteria Met)
Criteria_Met_LP_ON_data_raw_F = raw_F_LP_ON(GCAMP.Criteria_met,:);
%z-score across trials on all lps, then segment them out by met vs fail
raw_F_LP_ON_Z = zscore(raw_F_LP_ON);
raw_F_LP_OFF_Z = zscore(raw_F_LP_OFF);
Criteria_Met_LP_ON_data_Z_raw = raw_F_LP_ON_Z(GCAMP.Criteria_met, :);
Criteria_Met_LP_OFF_data_Z_raw = raw_F_LP_OFF_Z(GCAMP.Criteria_met, :);
Criteria_Fail_LP_ON_data_Z_raw = raw_F_LP_ON_Z(GCAMP.Criteria_fail, :);
Criteria_Failt_LP_OFF_data_Z_raw = raw_F_LP_OFF_Z(GCAMP.Criteria_fail, :);

%% Non-Reinforced Lever Presses (Criteria Fail)
Criteria_Fail_LP_ON_data_raw_F = raw_F_LP_ON(GCAMP.Criteria_fail, :);
%% Plot the average of criteria (met vs fail) lever press onset
plot_time = base_time_end:1/GCAMP.SR:time_end;

%% Raw f criteria
y7 = Criteria_Met_LP_ON_data_raw_F;
y8 = Criteria_Fail_LP_ON_data_raw_F;
n7 = size(Criteria_Met_LP_ON_data_raw_F,1);
n8 = size(Criteria_Fail_LP_ON_data_raw_F,1);
y9 = raw_F_LP_ON;
n9 = size(raw_F_LP_ON,1);

figure('Name',[GCAMP.mouseID ' ' GCAMP.training_day ' Criteria Raw F'],'NumberTitle','off', 'rend','painters','pos',[10 10 1200 700]);
%figure('rend','painters','pos',[10 10 1200 700])
subplot(2,2,[1 3])
hold on

%Plot mean with shaded standard error 
s = shadedErrorBar(plot_time, y7, {@mean, @(x) std(x) / sqrt(n7)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

s = shadedErrorBar(plot_time, y8, {@mean, @(x) std(x) / sqrt(n8)}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

% Plot shaded criteria time
zl = ylim;
% Create custom legend (to ignore shaded regions)
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-b');
h(2) = plot(NaN,NaN,'-r');
legend(h,{'Criteria Met', 'Criteria Fail'})
set(h,'LineWidth',4);
legend boxoff
xlabel('Time from Lever Press Onset (S)')
ylabel('Raw Intensity', 'FontWeight','bold')
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')
hold off

%% Heatmap
[out,idx] = sort(GCAMP.HoldDown_times);
data = raw_F_LP_ON(idx,:);

subplot(2,2,[2 4])
minValue = min(data(:));
maxValue = max(data(:));
% Scale the data to between -1 and +1.
data = (data-minValue) * 2 / (maxValue - minValue) - 1;

imagesc(data);
colormap('viridis')
h = colorbar;
ylabel(h,'Normalized \Delta F / F (%)', 'FontWeight','bold')
xticks(1:GCAMP.SR:length(plot_time));

title('\Delta F / F (%) for each Lever Press')
xlabel('Time from Lever Press Onset (S)')
ylabel('Lever Presses (Descending in Length)', 'FontWeight','bold')
xticklabels({'-5', '-4', '-3', '-2','-1', '0', '1','2','3','4', '5'})
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')

% Overlay Onset, Criteria, and LP Offsets
adjustedtimes = GCAMP.HoldDown_times(idx)/((1/GCAMP.SR)*1000);
shifted = adjustedtimes + (GCAMP.SR*abs(base_time_end)+1);
hold on
plot(shifted,1:length(idx), '-r','LineWidth', 2);
shiftedCriteria = (GCAMP.Criteria/((1/GCAMP.SR)*1000)) + (GCAMP.SR*abs(base_time_end)+1);
shiftedCriteriaMatrix = 1:length(idx);
shiftedCriteriaMatrix(:,:) = shiftedCriteria;
plot(shiftedCriteriaMatrix,1:length(idx), 'w','LineWidth',2);
zeroline = 1:length(idx); 
zeroline(:,:) = (GCAMP.SR*abs(base_time_end)+1);
plot(zeroline,1:length(idx),'w','LineWidth',2);
hold off

%% Presses in order of occurrence, raw signal
figure('Name',[GCAMP.mouseID ' ' GCAMP.training_day ' Criteria Raw F in serial order'],'NumberTitle','off', 'rend','painters','pos',[10 10 1200 700]);

[out,idx] = sort(GCAMP.HoldDown_times);
data = raw_F_LP_ON;

imagesc(data);
colormap('viridis')
h = colorbar;
ylabel(h,'Raw Intensity', 'FontWeight','bold')
xticks(1:GCAMP.SR:length(plot_time));

title('Raw Intensity,Order of occurrence')
xlabel('Time from Lever Press Onset (S)')
ylabel('Lever Presses(serial order)', 'FontWeight','bold')
xticklabels({'-5', '-4', '-3', '-2','-1', '0', '1','2','3','4', '5'})
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')

% Overlay Onset, Criteria, and LP Offsets
adjustedtimes = GCAMP.HoldDown_times/((1/GCAMP.SR)*1000);

shifted = adjustedtimes + (GCAMP.SR*abs(base_time_end)+1);
hold on
plot(shifted,1:length(idx), '-r','LineWidth', 2);
shiftedCriteria = (GCAMP.Criteria/((1/GCAMP.SR)*1000)) + (GCAMP.SR*abs(base_time_end)+1);
shiftedCriteriaMatrix = 1:length(idx);
shiftedCriteriaMatrix(:,:) = shiftedCriteria;
plot(shiftedCriteriaMatrix,1:length(idx), 'w','LineWidth',2);
zeroline = 1:length(idx); 
zeroline(:,:) = (GCAMP.SR*abs(base_time_end)+1);
plot(zeroline,1:length(idx),'w','LineWidth',2);
hold off


%% raw zscore baseline graphs
figure('Name',[GCAMP.mouseID ' ' GCAMP.training_day 'Raw-Z-score to baseline'],'NumberTitle','off','rend','painters','pos',[10 10 1200 700])
subplot([121])
hold on
s = shadedErrorBar(plot_time, smoothdata(raw_baseline_z_score_Met,2,'gaussian',5), {@mean, @(x) std(x) / sqrt(size(raw_baseline_z_score_Met,1))}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];
% 
s = shadedErrorBar(plot_time, smoothdata(raw_baseline_z_score_Fail,2,'gaussian',5), {@mean, @(x) std(x) / sqrt(size(raw_baseline_z_score_Fail,1))}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

zl = ylim;
CritInSecs = GCAMP.Criteria/1000;
ha = area([0 CritInSecs], [zl(1) zl(1)], 'FaceAlpha', 0.30, 'EdgeAlpha', 0);
he = area([0 CritInSecs], [zl(2) zl(2)], 'FaceAlpha', 0.30, 'EdgeAlpha', 0);
title(['Z-Scored to Baseline (raw)']);
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-b');
h(2) = plot(NaN,NaN,'-r');
set(h,'LineWidth',4);
legend(h,{'Met', 'Fail'},'FontSize', 20)
legend boxoff
xlabel('Time from Lever Press Onset (S)')
ylabel('Z-Score', 'FontWeight','bold')
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')
hold off

%% Heatmap

[out,idx] = sort(GCAMP.HoldDown_times);
data = baseline_z_score_LP_On_Raw(idx,:);

subplot(2,2,[2 4])
minValue = min(data(:));
maxValue = max(data(:));
% Scale the data to between -1 and +1.
data = (data-minValue) * 2 / (maxValue - minValue) - 1;

imagesc(data);
colormap('viridis')
h = colorbar;
ylabel(h,'Normalized \Delta F / F (%)', 'FontWeight','bold')
xticks(1:GCAMP.SR:length(plot_time));

title('\Delta F / F (%) for each Lever Press')
xlabel('Time from Lever Press Onset (S)')
ylabel('Lever Presses (Descending in Length)', 'FontWeight','bold')
xticklabels({'-5', '-4', '-3', '-2','-1', '0', '1','2','3','4', '5'})
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')

% Overlay Onset, Criteria, and LP Offsets
adjustedtimes = GCAMP.HoldDown_times(idx)/((1/GCAMP.SR)*1000);
shifted = adjustedtimes + (GCAMP.SR*abs(base_time_end)+1);
hold on
plot(shifted,1:length(idx), '-r','LineWidth', 2);
shiftedCriteria = (GCAMP.Criteria/((1/GCAMP.SR)*1000)) + (GCAMP.SR*abs(base_time_end)+1);
shiftedCriteriaMatrix = 1:length(idx);
shiftedCriteriaMatrix(:,:) = shiftedCriteria;
plot(shiftedCriteriaMatrix,1:length(idx), 'w','LineWidth',2);
zeroline = 1:length(idx); 
zeroline(:,:) = (GCAMP.SR*abs(base_time_end)+1);
plot(zeroline,1:length(idx),'w','LineWidth',2);
hold off



%% Interpolation
interp_size = 20; % how many samples you want to map activity to
interp_data = nan(length(GCAMP.HoldDown_times), interp_size); % placeholder matrix for interpolated data
zero_index = 101;
max_length = zero_index + time_end*GCAMP.SR;
timeinsamples = round((GCAMP.HoldDown_times)/((1/GCAMP.SR)*1000));
for num_trial = 1:length( GCAMP.HoldDown_times) % for each lever press
        lever_press_offset_index = zero_index + timeinsamples(num_trial);
        
        if lever_press_offset_index > max_length
           lever_press_offset_index = max_length;
        end
        
        if  lever_press_offset_index >= zero_index +2
            
        interp_trial = baseline_z_score_LP_On_Raw(num_trial, zero_index:lever_press_offset_index); % variable vector of F (determined by duration) 
        total_samples = length(interp_trial);
        t0 = linspace(1,total_samples,total_samples); %original time vector
        t1 = linspace(1,total_samples,interp_size); % new time vector (specifying the time points at which you want to interpolate)
        interp_data(num_trial,:) = interp1(t0,interp_trial,t1,'makima'); %y0 data interpolated to fixed sample space
         end

end
% prior to plotting, remove any NaNs (which are LPs withouth at least 2
% samples)
graph_interp = rmmissing(interp_data);

percent_of_press = 5*(1:20);
%divide the interpolated  data into success and failed presses
%based on the 
interp_met =interp_data(GCAMP.Criteria_met,:);
interp_fail = interp_data(GCAMP.Criteria_fail,:);

graph_interp_met = rmmissing(interp_met);
graph_interp_fail = rmmissing(interp_fail);
n1 = size(graph_interp,1);
n2 = size(graph_interp_met,1);
n3 = size(graph_interp_fail,1);

figure('Name',[GCAMP.mouseID ' ' GCAMP.training_day ' Criteria Delta F'],'NumberTitle','off', 'rend','painters','pos',[10 10 1200 700]);
hold on

%Plot mean with shaded standard error 
s = shadedErrorBar(percent_of_press, graph_interp, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-k', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];
hold on

s = shadedErrorBar(percent_of_press, graph_interp_met, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];
hold on

s = shadedErrorBar(percent_of_press, graph_interp_fail, {@mean, @(x) std(x) / sqrt(n3)}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];
hold on

title(['Interpolated Press Durations'])

% Create custom legend (to ignore shaded regions)
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-k');
h(2) = plot(NaN,NaN,'-b');
h(3) = plot(NaN,NaN,'-r');
legend(h,{'All', 'Criteria Met', 'Criteria Fail'})
set(h,'LineWidth',4);
legend boxoff
xlabel('Percent of Press Duration')
ylabel('Z scored to baseline', 'FontWeight','bold')
set(gca,'FontSize',20)
set(gca, 'FontName', 'Arial')

GCAMP.graph_interp_met = graph_interp_met;
GCAMP.graph_interp_fail = graph_interp_fail;
GCAMP.interp_data = interp_data;
GCAMP.graph_interp = graph_interp;
 
%% Save data
criteria_percent = 100*(length(GCAMP.Criteria_met)/length(GCAMP.LP_ON_timestamps));
GCAMP.criteria_percent = criteria_percent;
GCAMP.base_time_start = base_time_start;
GCAMP.base_time_start = base_time_end;
GCAMP.time_end = time_end;
GCAMP.Closest_LP_ON_dist = Closest_LP_ON_dist_all;

GCAMP.Raw_LP_met = Criteria_Met_LP_ON_data_raw_F;
GCAMP.Raw_LP_fail =  Criteria_Fail_LP_ON_data_raw_F;
GCAMP.raw_baseline_z_score_Met = raw_baseline_z_score_Met;
GCAMP.raw_baseline_z_score_Fail = raw_baseline_z_score_Fail; 
GCAMP.baseline_z_score_LP_On_Raw = baseline_z_score_LP_On_Raw;
GCAMP.raw_F_LP_ON = raw_F_LP_ON;
GCAMP.raw_F_LP_OFF = raw_F_LP_OFF;
GCAMP.baseline_z_score_LP_OFF_Raw = baseline_z_score_LP_OFF_Raw;
GCAMP.baseline_z_score_LP_OFF_Raw_Met = baseline_z_score_LP_OFF_Raw_Met;
GCAMP.baseline_z_score_LP_OFF_Raw_Fail = baseline_z_score_LP_OFF_Raw_Fail;
GCAMP.baseline_z_score_LP_OFF_Raw = baseline_z_score_LP_OFF_Raw;
GCAMP.baselines_on_raw = All_LP_ON_data_Baselines_raw;
GCAMP.baselines_on_raw_met = All_LP_ON_data_Baselines_raw(GCAMP.Criteria_met,:);
GCAMP.baselines_on_raw_fail = All_LP_ON_data_Baselines_raw(GCAMP.Criteria_fail,:);
GCAMP.raw_F_LP_ON_Z =raw_F_LP_ON_Z;
GCAMP.raw_F_LP_OFF_Z = raw_F_LP_OFF_Z;
GCAMP.Criteria_Met_LP_ON_data_Z_raw = Criteria_Met_LP_ON_data_Z_raw;
GCAMP.Criteria_Met_LP_OFF_data_Z_raw =Criteria_Met_LP_OFF_data_Z_raw;
GCAMP.Criteria_Fail_LP_ON_data_Z_raw = Criteria_Fail_LP_ON_data_Z_raw;
GCAMP.Criteria_Failt_LP_OFF_data_Z_raw = Criteria_Failt_LP_OFF_data_Z_raw;
end

