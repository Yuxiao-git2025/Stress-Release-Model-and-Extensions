% =========================================================================
% Simulation with ASR model
% Modified 2026/03/19
% You can simulate one catalog, and then calculate its conditional
% intensity by function:<SRM_Estimate> which based on the SRM-Model
% =========================================================================
global Mmin Mmax Bvalue Tlim A B C eta
% Lenght of the catalog (years)
Tlim = 2^10 ;
% Minimum magnitude of the catalog (Mw)
Mmin = 4 ;
% Maximum energy stored in the fault or upper limit
Mmax = 7.5 ;
% B-value
Bvalue = 1.0 ;
% rate parameters
A=-1;
B=0.01;
C=0.5;
% Release efficiency
eta=0.75; 
clc;
[EventTime,EventMag,E,EventEnergy,OutSumEnergy]=ASR_SimCatalog_Barani(2025);
% Results=ASR_SimCatalog_Jaume(2025);

% =========================================================================
% plot results
ASR_PlotResults(EventTime,EventMag)