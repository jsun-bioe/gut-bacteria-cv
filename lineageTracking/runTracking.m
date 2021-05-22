% load Morphometrics output
s = load(fileName,'frame');
frame = s.frame;
load('colormap.mat')

% basic config
minArea = 100; %cell area threshold
linkDepth = 4; %maxmimum length in time to link 2 cells
cmArea = 0.3; %area overlap threshold between frames

%% drift correction
%frame = fixDrift(frame);

%% link cells up
% return with new fields 'check', 'framePointer', 'objectPointer' in 'object'
frame2 = createTree(frame, minArea,linkDepth, cmArea);

%% rearrange data into a tree structure
% return a trajectory table 'cells' and new fields 'count' and 'ID' in 'object'
[cells, frame3] = getTree(frame2);

%% visualize tree
% plot cell trajectories
visualizeTree(cells,0);

%% fix merging artifacts and visualize again
% config
fixDepth = 5;
fixArea = 0.2;

% fix missing and merging
cells2 = fixMerge(cells,frame3,fixDepth,fixArea);
visualizeTree(cells2,1);

%% save results
if ~exist(outputName)
    save(outputName,'cells','cells2','frame3','fixDepth','fixArea')
else
    fprintf(['File exist!','\n'])
end