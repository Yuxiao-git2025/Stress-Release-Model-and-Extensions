%% Script to obtain the intensity of SRM-Model
global evalNums
evalNums=200;
% [time,mag,Mmin]=SRM_GetSynthetic();
time=EventTime;
mag=EventMag;
Mmin=4;
% Compile the parameter A, B, C
A = -1.0; 
B = 0.01; 
C = 0.5;
% Calculate the conditional intensity
lambda=SRM_Estimate(time,mag,Mmin,[A,B,C]);
% AIC0=SRM_Estimate(time,mag,Mmin,[A,B,C],true);
SRM_PlotResults(lambda,time,mag,Mmin);

%% MLE-Fitting
[pNew,AIC]=SRM_MLE(time,mag,Mmin,[A,B,C]);
SRM_PlotResults(SRM_Estimate(time,mag,Mmin,pNew),time,mag,Mmin);

%% RPP-analysis
CumInt=SRM_Integral(time,mag,Mmin,pNew);
SRM_RPP(CumInt,time);
SRM_RPP_Transformed(CumInt,time);
