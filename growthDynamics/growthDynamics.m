%==== 2020/05/13, by Jiawei ====%
% extracting growth parameters  %
%===============================%
clear; clc; close all;

[FileName,PathName] = uigetfile('*_tracking.mat','Select the data file','MultiSelect','on');
if ~iscell(FileName)
    FileName = {FileName};
end

FIX_MERGE_ = 1;
fixDepth_default = 12;
pixel2micron = 0.064;
fixArea_default = 0.4;
dT = 1; % min

fData = {};
cData = {};
fN = [];
for ii=1:length(FileName)
    fileName = strcat(PathName,FileName{ii});
    outputName = strcat(fileName(1:end-4),'_dynamics.mat');
    fprintf(['Computing dynamics for time-lapse ',num2str(ii),' of ',...
        num2str(length(FileName)),' ...','\n']);
    computeDynamics
    fData{ii} = fdata;
    cData{ii} = cdata;
    fN = [fN length(fdata)];
end
fN = max(fN);
clear fdata cdata

%% plot dynamics
plotDynamics