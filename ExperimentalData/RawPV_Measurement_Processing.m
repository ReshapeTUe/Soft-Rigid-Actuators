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
measurementset = 'PV';
actuatortype = 'RBIE';

PlotFolder = 'Plots';

addpath([measurements '\' measurementset '\' actuatortype])

dinfo = dir([measurements '\' measurementset '\' actuatortype '\' '*.txt']);

filenames = {dinfo.name};

MeasureData = cell(length(filenames), 4);

%% Preliminary information internal volumes (all in mL)


V_Syringe = 20;      % Starting volume syringe (set at 20 mL)


%%%%%% Tubing averaged per actuator DETERMINE%%%%%%%     
V_tubes = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_FS = 6.02396+V_tubes;     % Internal volume fully soft actuator
    
V_intb = 0.11155;   % Internal shells, big
V_ints = 0.05488;   % Internal shells, smell
       


V_RBIEI = 1.270251747;  % Internal volume tiny actuator segment (no inlet)

V_IL = 0.057176986;     % Long inlet internal volume
V_IS = 0.035342917;     % Short inlet internal volume
V_IAS = 0.007952156;    % Inter-Actuator tubing space



V_FB = V_FS;
V_FBI = V_FS - 13*V_intb - 2*V_ints;
V_FSE = V_FS;
V_FBIE = V_FS - 13*V_intb - 2*V_ints;
V_RBIE = 3*V_RBIEI + V_IL + 4*V_IS + 2*V_IAS - 9*V_intb;



%% Continue from here -->



%% Processing

for i=1:length(filenames)
   TableData = readtable(filenames{i});
   MatrixData = TableData{:,:};
   MatrixData = MatrixData - MatrixData(1,:);
       
   MeasureData{i,1} = filenames{i}(1:end-4);
   MeasureData{i,2} = MatrixData(:,2);      % Pressure
   MeasureData{i,4} = MatrixData(:,1);      % Time
  
   %%%% Add Volume Conversion Here
   MeasureData{i,3} = MatrixData(:,3);      % Reduced Syringe Volume
   
   MaxV(i) = max(MatrixData(:,3));
   
   clear MatrixData TableData
end 

%% Plotting 

Fig = figure();
hold on
grid on
box on
for i = 1:length(filenames)
    plot(MeasureData{i,3}, MeasureData{i,2}, 'Linewidth', lw)
    maxVol(i) = max(MeasureData{i,3});
    maxPressure(i)= max(MeasureData{i,2});
    V_min(i) = min(MeasureData{i,3});
    
end
xlabel('Pressure (KPa)')
ylabel('Force (N)')
if min(V_min) < -10
     xlim([min(V_min) 0])
     ylim([0 max(maxPressure)+0.05])
     legend(MeasureData{:,1}, 'Location', 'northeast')
else
    xlim([0 max(maxVol)])
    ylim([0 max(maxPressure)+0.3])
    legend(MeasureData{:,1}, 'Location', 'northwest')
end



