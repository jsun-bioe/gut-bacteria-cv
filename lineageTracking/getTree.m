function [cells, frame] = getTree(frame)
% Organize data into a tree structure using the pointers

% create field 'count' and 'ID'
for i=1:length(frame)
    for j=1:length(frame(i).object)
        frame(i).object(j).count = 0; %Counting the number of times getting pointed at (number of descendants)
        frame(i).object(j).ID = 0; %ID for cell trajectories in 'cells'
    end
end

% every cell gives a count to its mother cell(s)
for i=length(frame):-1:2 %loop backwards over frames
    object = frame(i).object;
    for j=1:length(object)
        if object(j).check==1 %passed area filter
            framePointer = object(j).framePointer;
            objectPointer = object(j).objectPointer;
            for k=1:length(framePointer)
                if framePointer(k)==4 && objectPointer(k)==2
                    pause
                end
                frame(framePointer(k)).object(objectPointer(k)).count =...
                frame(framePointer(k)).object(objectPointer(k)).count + 1;
            end
        end
    end
end

% stitch together cell trajectories with the counts
% 'cells' contains fields: frame, object, mother, daughter, merge
cells = [];
for i=length(frame):-1:1%loop backwards over frames
    object = frame(i).object;
    
    % last frame only
    if i==length(frame)
        for j=1:frame(i).num_objs
            if object(j).check==1
                cell = [];
                cell.frame = i;
                cell.object = j;
                cell.mother = []; %empty: unknown or does not exist
                cell.daughter = []; %empty: does not exist
                cell.merge = 0; %merge: ends by merging to another cell; 0: does not merge
                cells = [cells cell];
                % assign a trajectory ID
                object(j).ID = length(cells); %0: does not appear in 'cells'
            end
        end
    end
    
    % trace back trajectories by looking for mother(s)
    for j=1:length(object)
        if object(j).check==1
            
            framePointer = object(j).framePointer;
            objectPointer = object(j).objectPointer;
            ID = object(j).ID;
            
            % initiate a trajectory
            if ID==0 %end of a trajectory
                cell = [];
                cell.frame = i;
                cell.object = j;
                cell.mother = [];
                cell.daughter = [];
                cell.merge = 0;
                cells = [cells cell];
                % assign ID
                ID = length(cells);
                object(j).ID = ID;
            end
            
            % growth the trajectory by adding mother(s)
            if length(framePointer)>1 %a merged cell (artifact)
                % merging/division is treated as an end of a trajectory
                for k=1:length(framePointer)
                    mother_ID = frame(framePointer(k)).object(objectPointer(k)).ID;
                    if mother_ID==0 %initiate a new trajectory (a new entry in 'cells')
                        cell = [];
                        cell.frame = framePointer(k);
                        cell.object = objectPointer(k);
                        cell.mother = [];
                        cell.daughter = [];
                        cell.merge = ID;
                        cells = [cells cell];
                        % assign ID
                        frame(framePointer(k)).object(objectPointer(k)).ID = length(cells);
                        % add mother
                        cells(ID).mother = [cells(ID).mother length(cells)];
                    else %if entry exist (maybe already in a trajectory, e.g. A->(C,D), B->D)
                        cells(mother_ID).merge = [cells(mother_ID).merge ID];
                        cells(ID).mother = [cells(ID).mother mother_ID];
                    end
                end
                clear mother_ID
            elseif length(framePointer)==1 % just 1 mother (normal)
                if frame(framePointer).object(objectPointer).count==1 %just previous self
                    cells(ID).frame = [framePointer cells(ID).frame];
                    cells(ID).object = [objectPointer cells(ID).object];
                    % assign ID
                    frame(framePointer).object(objectPointer).ID = ID;
                elseif frame(framePointer).object(objectPointer).count>1 %real mother cell
                    mother_ID = frame(framePointer).object(objectPointer).ID;
                    if mother_ID==0 %create new entry (division is treated as an end of trajectory)
                        cell = [];
                        cell.frame = framePointer;
                        cell.object = objectPointer;
                        cell.mother = [];
                        cell.daughter = ID;
                        cell.merge = [];
                        cells = [cells cell];
                        % assign ID
                        frame(framePointer).object(objectPointer).ID = length(cells);
                        % add mother
                        % NOTE that trajectories will not have unempty
                        %'merge' and 'daughter' fields at the same time!
                        cells(ID).mother = length(cells);
                    else %entry exist
                        cells(mother_ID).daughter = [cells(mother_ID).daughter ID];
                        cells(ID).mother = mother_ID;
                    end
                    clear mother_ID
                end
            elseif isempty(framePointer) %starter (can be artifact)
                % no action to take
            end
        end
    end
    frame(i).object = object;
end

end

