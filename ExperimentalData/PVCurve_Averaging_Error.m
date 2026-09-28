clear all; clc
close all;
%% Set text interperters to 'latex'
set(0, 'defaultTextInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Set Figure Format
textwidth = 7;
textheight = 7;
fsize = 12; % Font size
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


V_max = 20;


V_step = 2;
Vrange = 0:sign(V_max)*V_step:V_max;
dV = V_step/2;

P_cur = NaN*zeros(length(filenames), length(Vrange));

for j = 1:length(filenames)
    for i = 1:length(Vrange)
        V_cur = Vrange(i);
        P_Raw = MeasureData{j,2};
        V_raw = MeasureData{j,3};
        
        P_window = P_Raw(V_raw>=V_cur-dV & V_raw<V_cur+dV);
        
        P_cur(j,i) = mean(P_window,"omitnan");
    end
end

P_cur(:,1) = 0;
P_mean = mean(P_cur, "omitnan");


% N* std = n confidence interval (1 = 66, 2 = 95, 3 = 99.7)
SD = 3*std(P_cur, "omitnan");


SE = zeros(length(SD),1);
for i = 1:length(Vrange)
    Nsamples = nnz(P_cur(:,i)) - nnz(isnan(P_cur(:,i)));
    SE(i) = SD(i)/sqrt(Nsamples);
end

Fig = figure(1);
hold on;
box on;
AVG = errorbar(Vrange, P_mean, SE, 'o', 'Linewidth', lw, 'MarkerFaceColor',[0.65 0.85 0.90]);
set(gca,'fontsize', fsize)
set(gca,'TickLength',[0.035 0.035]);
xlabel("\textbf{Reduced Syringe Volume ($\mathbf{mL}$)}", 'FontSize', fsize)
ylabel("\textbf{Pressure ($\mathbf{KPa}$)}", 'FontSize', fsize)

xlim([-2 22])
ylim([-5 90])
yticks([0 20 40 60 80 100])
xticks(0:5:20)

max(P_mean)










