%=== 2020/04/29, by Jiawei ===%
%    cell lineage tracking    %
%=============================%
clear; clc; close all

[FileName,PathName] = uigetfile('*_CONTOURS.mat','Select the data file','MultiSelect','on');

if ~iscell(FileName)
    FileName = {FileName};
end
for i=1:length(FileName)
    fileName = strcat(PathName,FileName{i});
    outputName = strcat(fileName(1:end-4),'_tracking.mat');
    fprintf(['Tracking time-lapse ',num2str(i),' of ',num2str(length(FileName)),' ...','\n']);
    runTracking
end