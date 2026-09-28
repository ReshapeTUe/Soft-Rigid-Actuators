clear all; close all; clc

%% Set text interperters to 'latex'
set(0, 'defaultTextInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Set Figure Format
textwidth = 7;
textheight = 7;
fsize = 12; % Font size
lw = 1.5;

%% Data

Separated = 360 + [239.23 258.35 217.75 216.95 272.86];
Bridge = [370.77 311.70 335.81 330.17 335.17 355.99];
Bridge_Internal = [300.75 313.96 295.00 310.46 321.50];
Outer = 360+ [236.12 226.49 241.65 216.44 259.88 265.29];
OuterInnerBridge = 360 + [105.68 125.62 148.69 125.64 127.62];
Skeleton = [183.07 179.47 205.12 176.81 175.93 181.56];


%% Processing



Set_Cur = NaN*zeros(6);

Set_Cur(1:5, 1) = Separated(1:5)';
Set_Cur(1:5, 2) = Bridge(1:5)';
Set_Cur(1:5, 3) = Bridge_Internal(1:5)';
Set_Cur(1:5, 4) = Outer(1:5)';
Set_Cur(1:5, 5) = OuterInnerBridge(1:5)';
Set_Cur(1:5, 6) = Skeleton(1:5)';


% N* std = n confidence interval (1 = 66, 2 = 95, 3 = 99.7)
Angle_mean = mean(Set_Cur, "omitnan");

SD = 3*std(Set_Cur, "omitnan");

SE = zeros(length(SD),1);

for i = 1:6
    Nsamples(i) = nnz(Set_Cur(:,i)) - nnz(isnan(Set_Cur(:,i)));
    SE(i) = SD(i)/sqrt(Nsamples(i));
end

TextLabels = categorical(["\textbf{(A)}", "\textbf{(B)}", "\textbf{(C)}", "\textbf{(D)}", "\textbf{(E)}", "\textbf{(F)}"]);




Fig = figure(1);
hold on;
grid on;
box on;
bar(TextLabels, Angle_mean)
xline(3.5, 'k--', 'Linewidth', 2*lw)
er = errorbar(TextLabels, Angle_mean, SE, 'Linewidth', lw);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

set(gca,'fontsize', fsize-3)
xlabel("Actuator Type", 'FontSize', fsize)
ylabel("Angular Displacement at Max Pressure (deg)", 'FontSize', fsize)








