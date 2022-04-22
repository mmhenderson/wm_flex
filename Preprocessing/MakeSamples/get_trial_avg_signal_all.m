% Helper script to compute trial-averaged signal in each voxel in each
% task, for desired subjects.
% these functions each load files named like "SampleFile_S0X.mat" and save
% files named like "MainTaskSignalByTrial_S0X.mat". 
% those files will all be in a directory called "Samples", which should be
% on same level as the directory "Preprocessing" (which is one level up
% from this directory).

clear
close all

sublist = [2:7];

getAvgSignal_mainTask(sublist)
getAvgSignal_SWMLoc(sublist)
getAvgSignal_spatLoc(sublist)
getAvgSignal_digLoc(sublist)
