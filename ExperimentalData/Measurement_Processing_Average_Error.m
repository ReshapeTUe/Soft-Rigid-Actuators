clear all;  clc
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
measurementset = 'Vacuum';
actuatortype = 'RBIE';

PlotFolder = 'Plots';

addpath([measurements '\' measurementset '\' actuatortype])

dinfo = dir([measurements '\' measurementset '\' actuatortype '\' '*.mat']);
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

if max(MeasureData{1,2})>70
   P_max = 100;
elseif min(MeasureData{1,2})< -20
   P_max = -80;
else
   P_max = 70;
end

P_step = 10;
Prange = 0:sign(P_max)*P_step:P_max;
dp = P_step/2;

F_cur = NaN*zeros(length(filenames), length(Prange));

for j = 1:length(filenames)
    
    for i = 1:length(Prange)
        P_cur = Prange(i);
        Force_Raw = MeasureData{j,4};
        P_raw = MeasureData{j,3};
        
        Force_window = Force_Raw(P_raw>=P_cur-dp & P_raw<P_cur+dp);
        
        F_cur(j,i) = mean(Force_window,"omitnan");
    end
end

F_cur(:,1) = 0;
F_mean = mean(F_cur, "omitnan");


% N* std = n confidence interval (1 = 66, 2 = 95, 3 = 99.7)
SD = 3*std(F_cur, "omitnan");


SE = zeros(length(SD),1);
for i = 1:length(Prange)
    Nsamples = nnz(F_cur(:,i)) - nnz(isnan(F_cur(:,i)));
    SE(i) = SD(i)/sqrt(Nsamples);
end

Fig = figure(1);
hold on;
box on;
if P_max == 70
     HatchX = [70 70 110 110];
     HatchY = [-0.2 2.4 2.4 -0.2];
     Env = polyshape(HatchX,HatchY);
     plot (Env, 'FaceColor',[0.7 0.7 0.7], 'FaceAlpha', 0.5)
end
AVG = errorbar(Prange, F_mean, SE, 'o', 'Linewidth', lw, 'MarkerFaceColor',[0.65 0.85 0.90]);
set(gca,'fontsize', fsize)
set(gca,'TickLength',[0.035 0.035]);
xlabel("\textbf{Pressure ($\mathbf{KPa}$)}", 'FontSize', fsize)
ylabel("\textbf{Force ($\mathbf{N}$)}", 'FontSize', fsize)
 







if P_max > 0
    xlim([-5 105])
    ylim([-0.1 2.4])
    xticks([0 20 40 60 80 100])
else
    xlim([-85 5])
    set ( gca, 'xdir', 'reverse' )
    ylim([-0.1 1.05])
    xticks([-80 -60 -40 -20 0])

end







