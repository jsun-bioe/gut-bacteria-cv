function frame = createTree(frame, minArea, linkDepth,cmArea)
% Cell tracking based on contour overlap

warningid = 'MATLAB:polyshape:repairedBySimplify';
warning('off', warningid);

% parameter (distance threshold in finding mother cell)
Dmax = 100; %unit pixel

% cell area filter
fprintf(['>> filtering...','\n'])
for i=1:length(frame)
    object = frame(i).object;
    for j=1:frame(i).num_objs
        if object(j).area < minArea
            object(j).check = 0; %0: did not pass area filter
        else
            object(j).check = 1; %clear for further analysis
        end
    end
    frame(i).object = object; %update object
end
clear I
fprintf(['>> filtering done!','\n'])

% cell tracking based on contour overlap
% framePointer: a field in 'object' pointing to the frame number of the cell
%at its last occurrence ("mother")
% objectPointer: a field in 'object' giving the specific object
%corresponding to this cell in the frame specified by framePointer
fprintf(['>> start tracking...','\n'])
for i=1:length(frame)
    % print progress
    if mod(i,10)==0
        fprintf(['   --- frame # ',num2str(i),' ---','\n'])
    end
    
    num_objs = frame(i).num_objs;
    object = frame(i).object;
    
    % frame #1
    if i==1
        for j=1:num_objs
            if object(j).check==1
                object(j).framePointer = []; %[]: no mother found yet
                object(j).objectPointer = [];
            else
                object(j).framePointer = -1; %-1: did not pass area filter
                object(j).objectPointer = -1;
            end
        end
        frame(1).object = object;
        continue
    end
    
    % frame # > 1
    frm = (max(1,i-linkDepth):1:i-1);
    for j=1:num_objs
        if object(j).check==0
            object(j).framePointer = -1;
            object(j).objectPointer = -1;
            continue
        end
        % find mother cell(s) for frame i, object j using area overlap
        cont1 = polyshape(object(j).Xcont, object(j).Ycont);
        framePointer = []; %empty: did not find mother (default)
        objectPointer = [];
        % record cells that have already been included in the tree of the current 
        %cell (cont1) by links created in previous loops to avoid overcounting
        pointed = [];
        for k=length(frm):-1:1 %trace back to find mother
            object2 = frame(frm(k)).object;
            for jj=1:frame(frm(k)).num_objs %loop over all contours
                if object2(jj).check==0 %did not pass area filter
                    continue
                end
                % filter based on distance (to reduce computing)
                if sqrt((object2(jj).Xcent-object(j).Xcent)^2+(object2(jj).Ycent-object(j).Ycent)^2) > Dmax
                    continue
                end
                cont2 = polyshape(object2(jj).Xcont, object2(jj).Ycont);
                overlap = intersect(cont1,cont2); %find overlap
                if (area(overlap)/area(cont1)>cmArea) || (area(overlap)/area(cont2)>cmArea)
                    if ~ismember(jj,pointed) %objects that are not included in the current tree
                        framePointer = [framePointer frm(k)];
                        objectPointer = [objectPointer jj];
                    end
                end
            end
            % IMPORTANT: trance back by 1 frame to find cells that are
            %already in current tree to avoid double counting
            if k>1
                pointed = traceBack(framePointer,objectPointer,frame,frm(k-1));
            end
        end
        % update mother cell information
        object(j).framePointer = framePointer;
        object(j).objectPointer = objectPointer;
    end
    frame(i).object = object;
end
fprintf(['>> tracking done!','\n'])

end

