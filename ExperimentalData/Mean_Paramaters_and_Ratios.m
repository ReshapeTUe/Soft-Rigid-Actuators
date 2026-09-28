close all;
clear all;
clc

%%%%%%%%%%%%%%% Means %%%%%%%%%%%%%%%%%

MeanFPos = [0.3307, 0.3253, 0.3791, 0.3934, 0.4883, 2.1780];
MeanFNeg = [0.0504, 0.1293, 0.3041, 0.0637, 0.3482, 0.7613];
MeanPV = [45.3524, 49.6518, 55.3757, 66.7866, 68.6948, 83.7149];
MeanAngle = [601.0280,336.724,308.334 ,596.116 ,486.65,184];

%%%%%%%%%%%%%% Percentages %%%%%%%%%%%%%%
PercFPos = 100*(MeanFPos/MeanFPos(1));
PercFNeg = 100*(MeanFNeg/MeanFNeg(1));
PercPV = 100*(MeanPV/MeanPV(1));
PercAngle = 100*(MeanAngle/MeanAngle(1));
