function [] = GCAMP_plot_with_baseline_Z_scores_r(Grouped_GCAMP, base_time_start, base_time_end,time_end)
%Last Updated 5/10/21 by Drew Schreienr
%takes the data from individual mice/sessions and puts it together
%as the data is baseline z-scored, we can have all of it on the same scale

Baseline_Z_Scored_Met =[];
Baseline_Z_Scored_Fail = [];
Baseline_Z_Scored_All =[];
z_interp_data =[];
z_interp_data_met =[];
z_interp_data_fail = [];
all_HD_times =[];
LP_Durations_All = [];
mean_duration =[];
Font_Size = 30;
smoothness = 10;
Baseline_Z_Scored_OFF_Met =[];
Baseline_Z_Scored_OFF_Fail = [];
blue = [0, 0.4470, 0.7410];
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
purple = [0.4940, 0.1840, 0.5560];
Colors = {blue, red, green, purple};

for i = 1:length(Grouped_GCAMP.Mice)
  
    Baseline_Z_Scored_Met =[Baseline_Z_Scored_Met Grouped_GCAMP.Mice{i}.raw_baseline_z_score_Met';];
    Baseline_Z_Scored_Fail = [Baseline_Z_Scored_Fail Grouped_GCAMP.Mice{i}.raw_baseline_z_score_Fail';];
    Baseline_Z_Scored_All = [Baseline_Z_Scored_All Grouped_GCAMP.Mice{i}.baseline_z_score_LP_On_Raw';];  
    Baseline_Z_Scored_OFF_Met =[Baseline_Z_Scored_OFF_Met Grouped_GCAMP.Mice{i}.baseline_z_score_LP_OFF_Raw_Met';];
    Baseline_Z_Scored_OFF_Fail = [Baseline_Z_Scored_OFF_Fail Grouped_GCAMP.Mice{i}.baseline_z_score_LP_OFF_Raw_Fail';];
    z_interp_data =[z_interp_data Grouped_GCAMP.Mice{i}.graph_interp';];
    z_interp_data_met =[z_interp_data_met Grouped_GCAMP.Mice{i}.graph_interp_met';];
    z_interp_data_fail = [z_interp_data_fail Grouped_GCAMP.Mice{i}.graph_interp_fail';];
       
    all_HD_times = [all_HD_times; Grouped_GCAMP.Mice{i}.HoldDown_times;];
    mean_duration = [mean_duration mean(Grouped_GCAMP.Mice{i}.HoldDown_times);];
    LP_Durations_All = [LP_Durations_All; Grouped_GCAMP.Mice{i}.HoldDown_times];
   
end

%re-arrange the individual aniaml z-score met v fail
%flip everything so its number of events by time
Baseline_Z_Scored_Met =Baseline_Z_Scored_Met';
Baseline_Z_Scored_Fail = Baseline_Z_Scored_Fail';
Baseline_Z_Scored_All = Baseline_Z_Scored_All';
Baseline_Z_Scored_OFF_Met =Baseline_Z_Scored_OFF_Met';
Baseline_Z_Scored_OFF_Fail = Baseline_Z_Scored_OFF_Fail';
z_interp_data =z_interp_data';
z_interp_data_met =z_interp_data_met';
z_interp_data_fail = z_interp_data_fail';

%% Baseline z-scored data

%lp onset
% CI of all lps
 bootCI = boot_CI(Baseline_Z_Scored_All,1000,.01);

%permutation test that requires 4 seq samples to pass significance
%threshold
 [p_val, observeddifference] = permTest_array(Baseline_Z_Scored_Met(:,41:133),Baseline_Z_Scored_Fail(:,41:133), 1000 );
 
threshold = .05;
    sample_span = 4;
aboveThreshold = p_val <= threshold;
        %aboveThreshold is a logical array, where 1 when above threshold, 0, below.
        %we thus want to calculate the difference between rising and falling edges
        aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
        edges = diff(aboveThreshold);
        rising = find(edges==1);     %rising/falling edges
        falling = find(edges==-1);
        spanWidth = falling - rising;  %width of span of 1's (above threshold)
        wideEnough = spanWidth >= sample_span;
        startPos = rising(wideEnough);    %start of each span
        endPos = falling(wideEnough)-1;   %end of each span
        %all points which are in the sample span (i.e. between startPos and endPos).
        allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));
        
   %now use that to index into the p_val thing and turn the significant
   %consecutive ones into 0s, and all others into Nans so that we can plot
   %that into the original graph to show where sig. diff exist.
   p_val_for_graph = p_val;
   p_val_for_graph(allInSpan) =0;
   non_sig_idx =p_val_for_graph > 0;
   p_val_for_graph(non_sig_idx) = NaN;
    
y1 = smoothdata(Baseline_Z_Scored_All,2,'gaussian',smoothness);
n1 = size(Baseline_Z_Scored_All ,1);
plot_time = base_time_end:1/Grouped_GCAMP.Mice{1}.SR:time_end;

mean_y1 = mean(Baseline_Z_Scored_All);
%need to subtrac t the mean from the CI, otherwise it plots mean + CI
bootCI_graph = [bootCI(1,:)-mean_y1;bootCI(1,:)-mean_y1] ;

figure('Name',[Grouped_GCAMP.Mice{1}.training_day 'Basline Z-scored All Lever Presses'],'NumberTitle','off','rend','painters','pos',[10 10 1200 850])
% subplot(2,2,[1 2])
hold on
%Plot mean with shaded standard error 
s = shadedErrorBar(plot_time,mean_y1, bootCI_graph, 'lineprops', '-k', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.5,0.25,0.25];
% 
zl = ylim;
ha = area([0 Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(1) zl(1)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);
he = area([0 Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(2) zl(2)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);

title(['Z-Scored Relative to Baseline'])
% Create custom legend (to ignore shaded regions)
h = zeros(1, 1);
h(1) = plot(NaN,NaN,'-k');
legend(h,{'All Lever Presses'})
set(h,'LineWidth',4);
legend boxoff
xlabel('Time from Lever Press Onset (Seconds)', 'FontName', 'Arial')
ylabel('Z scored Baseline', 'FontName', 'Arial')
set(gca,'FontSize',Font_Size)
set(gca, 'FontName', 'Arial')

y1 = smoothdata(Baseline_Z_Scored_Met(:,41:133),2,'gaussian',smoothness);
y2 = smoothdata(Baseline_Z_Scored_Fail(:,41:133),2,'gaussian',smoothness);
n1 = size(Baseline_Z_Scored_Met ,1);
n2 = size(Baseline_Z_Scored_Fail,1);

plot_time = base_time_end:1/Grouped_GCAMP.Mice{1}.SR:time_end;
plot_time = plot_time(:,41:133);

figure('Name',[Grouped_GCAMP.Mice{1}.training_day 'Basline Z-scored Criteria'],'NumberTitle','off','rend','painters','pos',[10 10 1200 850])
% subplot(2,2,[1 2])
hold on
%Plot mean with shaded standard error 
s = shadedErrorBar(plot_time, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.5,0.25,0.25];

s = shadedErrorBar(plot_time, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

% Plot shaded criteria time
zl = ylim;
ha = area([0 Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(1) zl(1)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);
he = area([0 Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(2) zl(2)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);

title(['Z-Scored Relative to Baseline'])
% Create custom legend (to ignore shaded regions)
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-b');
h(2) = plot(NaN,NaN,'-r');
legend(h,{'Met', 'Fail'})
set(h,'LineWidth',4);
legend boxoff
xlabel('Time from Lever Press Onset (Seconds)', 'FontName', 'Arial')
ylabel('Z scored Baseline',  'FontName', 'Arial')
set(gca,'FontSize',Font_Size)
set(gca, 'FontName', 'Arial')
xlim([-3 1.6])
 hold on 
 %this plots where the two different traces significantly differ at y = 0
 plot(plot_time,p_val_for_graph, 'k','LineWidth',4)
x0=10;
y0=10;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])

%Get mean and SD for export to prism
base_z_score_onset_met_mean =mean(y1,1);
base_z_score_onset_met_std = std(y1,1);
base_z_score_onset_fail_mean =mean(y2,1);
base_z_score_onset_fail_std = std(y2,1);

%% Baseline z-score for lp offset
 [p_val, observeddifference] = permTest_array(Baseline_Z_Scored_OFF_Met(:,69:end),Baseline_Z_Scored_OFF_Fail(:,69:end), 1000 );
 
threshold = .05;
    sample_span = 4;
aboveThreshold = p_val <= threshold;
        %aboveThreshold is a logical array, where 1 when above threshold, 0, below.
        %we thus want to calculate the difference between rising and falling edges
        aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
        edges = diff(aboveThreshold);
        rising = find(edges==1);     %rising/falling edges
        falling = find(edges==-1);
        spanWidth = falling - rising;  %width of span of 1's (above threshold)
        wideEnough = spanWidth >= sample_span;
        startPos = rising(wideEnough);    %start of each span
        endPos = falling(wideEnough)-1;   %end of each span
        %all points which are in the sample span (i.e. between startPos and endPos).
        allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));
        
   %now use that to index into the p_val thing and turn the significant
   %consecutive ones into 0s, and all others into Nans so that we can plot
   %that into the original graph to show where sig. diff exist.
   p_val_for_graph = p_val;
   p_val_for_graph(allInSpan) =0;
   non_sig_idx =p_val_for_graph > 0;
   p_val_for_graph(non_sig_idx) = NaN;
    
y1 = smoothdata(Baseline_Z_Scored_OFF_Met(:,69:end),2,'gaussian',smoothness);
y2 = smoothdata(Baseline_Z_Scored_OFF_Fail(:,69:end),2,'gaussian',smoothness);
n1 = size(Baseline_Z_Scored_OFF_Met ,1);
n2 = size(Baseline_Z_Scored_OFF_Fail,1);
plot_time = base_time_end:1/Grouped_GCAMP.Mice{1}.SR:time_end;
plot_time =plot_time(:,69:end);

%get mean and SD for prism graphing
base_z_score_offset_met_mean =mean(y1,1);
base_z_score_offset_met_std = std(y1,1);
base_z_score_offset_fail_mean =mean(y2,1);
base_z_score_offset_fail_std = std(y2,1);


figure('Name',[Grouped_GCAMP.Mice{1}.training_day 'Baseline Z-scored Criteria Offset'],'NumberTitle','off','rend','painters','pos',[10 10 1200 850])
% subplot(2,2,[1 2])
hold on
%Plot mean with shaded standard error 
s = shadedErrorBar(plot_time, y1, {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.5,0.25,0.25];

s = shadedErrorBar(plot_time, y2, {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

% Plot shaded criteria time
zl = ylim;
ha = area([0 -Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(1) zl(1)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);
he = area([0 -Grouped_GCAMP.Mice{1}.Criteria/1000], [zl(2) zl(2)], 'FaceAlpha', 0.20, 'EdgeAlpha', 0);
title(['Z-Scored to Baseline Press Offset'])
% Create custom legend (to ignore shaded regions)
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-b');
h(2) = plot(NaN,NaN,'-r');
legend(h,{'Criteria Met', 'Criteria Fail'})
set(h,'LineWidth',4);
legend boxoff
xlim([-1.6 5])
xlabel('Time from Lever Press Offset (Seconds)', 'FontName', 'Arial')
ylabel('Z scored Baseline', 'FontName', 'Arial')
set(gca,'FontSize',Font_Size)
set(gca, 'FontName', 'Arial')
hold on
 plot(plot_time,p_val_for_graph, 'k','LineWidth',4)

%% interpolated plots
[p_val, observeddifference] = permTest_array(z_interp_data_met,z_interp_data_fail, 1000 );
 
threshold = .05;
    sample_span = 3;
aboveThreshold = p_val <= threshold;
        %aboveThreshold is a logical array, where 1 when above threshold, 0, below.
        %we thus want to calculate the difference between rising and falling edges
        aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
        edges = diff(aboveThreshold);
        rising = find(edges==1);     %rising/falling edges
        falling = find(edges==-1);
        spanWidth = falling - rising;  %width of span of 1's (above threshold)
        wideEnough = spanWidth >= sample_span;
        startPos = rising(wideEnough);    %start of each span
        endPos = falling(wideEnough)-1;   %end of each span
        %all points which are in the sample span (i.e. between startPos and endPos).
        allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));
        
   %now use that to index into the p_val thing and turn the significant
   %consecutive ones into 0s, and all others into Nans so that we can plot
   %that into the original graph to show where sig. diff exist.
   p_val_for_graph = p_val;
   p_val_for_graph(allInSpan) =0;
   non_sig_idx =p_val_for_graph > 0;
   p_val_for_graph(non_sig_idx) = NaN;
    
n1 = size(z_interp_data ,1);
n2 = size(z_interp_data_met,1);
n3 = size(z_interp_data_fail,1);
percent_of_press = 5*(1:20);

figure('Name',[Grouped_GCAMP.Mice{1}.training_day 'Interpolated Press Durations'],'NumberTitle','off','rend','painters','pos',[10 10 1200 850])
hold on
%Plot mean with shaded standard error 
s = shadedErrorBar(percent_of_press, smoothdata(z_interp_data,2,'gaussian',5), {@mean, @(x) std(x) / sqrt(n1)}, 'lineprops', '-k', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.5,0.25,0.25];

s = shadedErrorBar(percent_of_press, smoothdata(z_interp_data_met,2,'gaussian',5), {@mean, @(x) std(x) / sqrt(n2)}, 'lineprops', '-b', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

s = shadedErrorBar(percent_of_press, smoothdata(z_interp_data_fail,2,'gaussian',5), {@mean, @(x) std(x) / sqrt(n3)}, 'lineprops', '-r', 'transparent',1);
set(s.edge,'LineWidth',1,'LineStyle','-')
s.mainLine.LineWidth = 4;
s.patch.FaceColor = [0.25,0.25,0.25];

title(['Interpolated Press Durations (z-scored to baseline)'])
% Create custom legend (to ignore shaded regions)
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-k');
h(2) = plot(NaN,NaN,'-b');
h(3) = plot(NaN,NaN,'-r');
legend(h,{'All', 'Criteria Met', 'Criteria Fail'})
set(h,'LineWidth',4);
legend boxoff
xlabel('Percent of Press Duration', 'FontName', 'Arial')
ylabel('Z scored to baseline', 'FontName', 'Arial')
set(gca,'FontSize',Font_Size)
set(gca, 'FontName', 'Arial')
hold on
plot(percent_of_press,p_val_for_graph, 'k','LineWidth',4)

 %mean and SD for prism graphing
base_interp_met_mean =mean( smoothdata(z_interp_data_met,2,'gaussian',5),1);
base_interp_met_std = std( smoothdata(z_interp_data_met,2,'gaussian',5),1);
base_interp_fail_mean =mean(smoothdata(z_interp_data_fail,2,'gaussian',5),1);
base_interp_fail_std = std(smoothdata(z_interp_data_fail,2,'gaussian',5),1);

end

