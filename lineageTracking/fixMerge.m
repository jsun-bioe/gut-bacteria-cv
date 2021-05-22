function cells = fixMerge(cells,frame,linkDepth,cmArea)
% Fix missing and merging artifacts and show information of trajectory connection with 
% a new field 'connected' in cells

warningid = 'MATLAB:polyshape:repairedBySimplify';
warning('off', warningid);
fprintf(['>> correcting artifacts...','\n'])

% parameter
searchLen = 100; %how far should we search in the 'cells' list for long-range missing

% clean field 'daughter' and 'merge' (get rid of the 0s)
for i=1:length(cells)
    cells(i).daughter = cells(i).daughter(cells(i).daughter~=0);
    cells(i).merge = cells(i).merge(cells(i).merge~=0);
end

% create a new field 'connetced'
% 0: not connected
% 1: connected to an upstream cell
% 2: connected to a downstream cell (if both, prioritize 1)
% -1: a merged cell being an artifact (corrected)
% -2: a merged cell being a true cell (no correction)
for i=1:length(cells)
    if length(cells(i).mother)>1
        cells(i).connected = -1; %default: artifact
    else
        cells(i).connected = 0; %default: no correction
    end
end

% fix long-range missing (larger than the previous linking depth)
for i=1:length(cells)
    if isempty(cells(i).mother) && cells(i).frame(1)>linkDepth
        cell0 = frame(cells(i).frame(1)).object(cells(i).object(1)); %first cell
        cont0 = polyshape(cell0.Xcont,cell0.Ycont);
        for j=i+1:1:min(length(cells),i+searchLen) %earlier trajectories (end point) appear later in cells
            if cells(j).frame(end) >= cells(i).frame(1) %earlier trajectory
                continue
            end
            if isempty(cells(j).daughter) %no daughter cell
                cell1 = frame(cells(j).frame(end)).object(cells(j).object(end));
                cont1 = polyshape(cell1.Xcont,cell1.Ycont);
                overlap = intersect(cont0,cont1); %also based on overlap
                if (area(overlap)/area(cont0)>cmArea) && (area(overlap)/area(cont1)>cmArea)
                    % trajectory i will then be linked to trajectory j
                    cells(j).frame = [cells(j).frame cells(i).frame];
                    cells(j).object = [cells(j).object cells(i).object];
                    cells(j).daughter = cells(i).daughter;
                    cells(j).merge = cells(i).merge;
                    cells(j).connected = 2;
                    cells(i).connected = 1;
                    % change the listed mothers of downstream daughter or
                    %merged cells as well
                    daughters = cells(i).daughter;
                    for k=1:length(daughters)
                        cells(daughters(k)).mother = j;
                    end
                    merging_to = cells(i).merge;
                    for k=1:length(merging_to)
                        cells(merging_to(k)).mother = j;
                    end
                    break %find at most one upstream cell (not looking for merging here)!
                end
            end
        end
    end
end
clear daughter merging_to

% fix merging artifact
for i=1:length(cells)
    if ((cells(i).connected==0) && ~isempty(cells(i).merge)) %&& (min(cells(i).merge)~=0)
        % get the ending cell
        cell0 = frame(cells(i).frame(end)).object(cells(i).object(end));
        % get merged starter cell
        merge_ID = cells(i).merge;
        
        % TO BE FIXED - merge_ID has more than 1 entry
        if length(merge_ID)>1
            merge_ID = merge_ID(1);
        end
        
        % the parameter 'linkDepth' sets a threshold for a trajectory being
        %real (long enough)
        if cells(merge_ID).frame(end)-cells(merge_ID).frame(1)+1>linkDepth
            % if a merged trajectory is considered not artifact, do nothing
            cells(merge_ID).connected = -2;
            continue
        end
        
        % if the merged cell is considered artifact
        cont0 = polyshape(cell0.Xcont,cell0.Ycont);
        % find daughter IDs
        daughters = cells(merge_ID).daughter;
        % check all daughters to see if they can be linked directly to
        %cell0
        for j=1:length(daughters)
            if (daughters(j)~=0)
                cell1 = frame(cells(daughters(j)).frame(1)).object(cells(daughters(j)).object(1));
                cont1 = polyshape(cell1.Xcont,cell1.Ycont);
                overlap = intersect(cont0,cont1);
                if (area(overlap)/area(cont0)>cmArea) && (area(overlap)/area(cont1)>cmArea) %notice this is &&
                    % connect frames & object
                    cells(i).frame = [cells(i).frame cells(daughters(j)).frame];
                    cells(i).object = [cells(i).object cells(daughters(j)).object];
                    % change labels
                    cells(i).connected = 2;
                    cells(daughters(j)).connected = 1;
                    % copy daughter & mother cells
                    cells(i).merge = cells(daughters(j)).merge;
                    for k=1:length(cells(i).merge)
                        cells(cells(i).merge(k)).mother = ...
                            [cells(cells(i).merge(k)).mother, i];
                        cells(cells(i).merge(k)).mother = ...
                            setdiff(cells(cells(i).merge(k)).mother, daughters(j));
                    end
                    cells(i).daughter = cells(daughters(j)).daughter;
                    for k=1:length(cells(i).daughter)
                        cells(cells(i).daughter(k)).mother = i;
                    end
                    
                    break % find at most one cell (may be improved in future)
                end
            end
        end
    end
end

% for 'merging mothers' that have no downstream trajectories, they are no
%longer a mother cell
for i=1:length(cells)
    if cells(i).connected == -2 % true merged cells
        mother = [];
        for j=1:length(cells(i).mother)
            if cells(cells(i).mother(j)).connected == 2
                mother = [mother j];
            end
        end
        cells(i).mother = mother;
    end
end

fprintf(['>> correction done!','\n'])

end