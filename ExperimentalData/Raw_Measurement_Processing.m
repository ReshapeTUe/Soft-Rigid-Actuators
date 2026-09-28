clear all; close all; clc

%% Set text interperters to 'latex'
set(0, 'defaultTextInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Set Figure Format
textwidth = 12;
textheight = 7;
fsize = 11; % Font size
lw = 1.5;

%% Read Data

measurements = 'measurements';
measurementset = 'Vacuum';
actuatortype = 'RBIE';

PlotFolder = 'Plots';

addpath([measurements '\' measurementtype '\' ActuatorType])

dinfo = dir([measurements '\' measurementtype '\' ActuatorType '\' '*.mat']);
filenames = {dinfo.name};

MeasureData = cell(length(filenames), 5);

%% Processing

for i=1:length(filenames)
load(filenames{i})
   p_des = p_des(~isnan(force_meas));
   p_meas = p_meas(~isnan(force_meas));
   time = time(~isnan(force_meas));
   force_meas = force_meas(~isnan(force_meas));

   
   for j = 2:length(force_meas)
       if abs(force_meas(j)-force_meas(j-1))>=0.5
           force_meas(j) = NaN;
       end
   end
   
   p_des = p_des(~isnan(force_meas));
   p_meas = p_meas(~isnan(force_meas));
   time = time(~isnan(force_meas));
   force_meas = force_meas(~isnan(force_meas));    
       
   MeasureData{i,1} = filenames{i}(1:end-4);
   MeasureData{i,2} = p_des;
   MeasureData{i,3} = p_meas;
   MeasureData{i,4} = force_meas-force_meas(1);
   MeasureData{i,5} = time-t_start;  
   

   clear p_des p_meas force_meas time t_start
end 

%% Plotting 

Fig = figure();
hold on
grid on
box on
for i = 1:length(filenames)
    plot(MeasureData{i,3}, MeasureData{i,4}, 'Linewidth', lw)
    maxpres(i) = max(MeasureData{i,3});
    maxforce(i)= max(MeasureData{i,4});
    
    P_min(i) = min(MeasureData{i,3});
    
end
xlabel('Pressure (KPa)')
ylabel('Force (N)')
if min(P_min) < -10
     xlim([min(P_min) 0])
     ylim([0 max(maxforce)+0.05])
     legend(MeasureData{:,1}, 'Location', 'northeast')
else
    xlim([0 max(maxpres)])
    ylim([0 max(maxforce)+0.3])
    legend(MeasureData{:,1}, 'Location', 'northwest')
end

