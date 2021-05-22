s = load(fileName);
cells = s.cells;
frame3 = s.frame3;
if isfield(s,'fixDepth')
    fixDepth = s.fixDepth;
    fixArea = s.fixArea;
else
    fixDepth = fixDepth_default;
    fixArea = fixArea_default;
end

%% fix merging
if FIX_MERGE_
    cells2 = fixMerge(cells,frame3,fixDepth,fixArea);
else
    cells2 = cells;
end

%% frame data
nframe = length(frame3);
fdata = [];
n0 = 3; % number of frames used to calculate growth rate
for i=1:nframe
    object = frame3(i).object;
    if isempty(object)
        continue
    end
    ID = cat(1,object.ID);
    
    % grwoth rate & area
    fdata(i).gr = [];
    fdata(i).A = [];
    fdata(i).dAdt = [];
    fdata(i).area = [];
    fdata(i).clen = [];
    fdata(i).SA = [];
    fdata(i).V = [];
    for j=1:length(ID)
        if FIX_MERGE_
            if (ID(j)~=0) && abs(cells2(ID(j)).connected)==1 % 1 or -1: should not be counted
                continue
            end
        end
        if (ID(j)~=0) && (object(j).check~=0)
            I0 = find(cells2(ID(j)).frame == i);
            I1 = find(cells2(ID(j)).frame == i+n0);
            if ~isempty(I1) % frame exists
                frm0 = i;
                frm1 = i+n0;
                obj0 = cells2(ID(j)).object(I0);
                obj1 = cells2(ID(j)).object(I1);
                a0 = frame3(frm0).object(obj0).area;
                a1 = frame3(frm1).object(obj1).area;
                fdata(i).gr = [fdata(i).gr (a1-a0)/n0/a0/dT]; % (dA/dt)/A
                fdata(i).dAdt = [fdata(i).dAdt (a1-a0)/n0/dT*pixel2micron^2]; % dA/dt
                fdata(i).A = [fdata(i).A a0*pixel2micron^2];
            end
            % area
            fdata(i).area = [fdata(i).area object(j).area*pixel2micron^2];
            % perimeter
            cont = polyshape(object(j).Xcont, object(j).Ycont);
            perim = perimeter(cont);
            fdata(i).clen = [fdata(i).clen perim*pixel2micron];
            frame3(i).object(j).perim = perim;
            % SA
            SA = surface_area_from_mesh(frame3(i).object(j).pill_mesh,...
                frame3(i).object(j).centerline);
            fdata(i).SA = [fdata(i).SA SA*pixel2micron^2];
            frame3(i).object(j).SA = SA;
            % V
            V = volume_from_mesh(frame3(i).object(j).pill_mesh,...
                frame3(i).object(j).centerline);
            fdata(i).V = [fdata(i).V V*pixel2micron^3];
            frame3(i).object(j).V = V;
        end
    end
    fdata(i).mean_area = mean(fdata(i).area);
    fdata(i).mean_clen = mean(fdata(i).clen);
end
clear ID I0 frm0 frm1 obj0 obj1 a0 a1 cont

%% cell data
ncell = length(cells);
cdata = [];
for i=1:ncell
    cdata(i).frame = cells2(i).frame;
    cdata(i).area = [];  % cell area
    cdata(i).clen = [];  % contour length
    cdata(i).SA = [];    % surface area
    cdata(i).V = [];     % volume
    if FIX_MERGE_ && (abs(cells2(i).connected)==1) % should not be counted
        cdata(i).bsize = []; % birth size
        cdata(i).dsize = []; % division size
        cdata(i).tb = [];    % time of birth
        cdata(i).td = [];    % time of division
    else
        if ~isempty(cells2(i).mother) && length(cells2(i).mother)==1 && cells2(i).mother>0
            cdata(i).bsize = frame3(cells2(i).frame(1)).object(cells2(i).object(1)).area*pixel2micron^2;
            cdata(i).tb = cells2(i).frame(1);
        elseif isempty(cells2(i).mother) && cells2(i).frame(1)<20 % starting cells
            cdata(i).bsize = frame3(cells2(i).frame(1)).object(cells2(i).object(1)).area*pixel2micron^2;
            cdata(i).tb = 1;
        else
            cdata(i).bsize = [];
            cdata(i).tb = [];
        end
        if ~isempty(cells2(i).daughter) && max(cells2(i).daughter)>0
            cdata(i).dsize = frame3(cells2(i).frame(end)).object(cells2(i).object(end)).area*pixel2micron^2;
            cdata(i).td = cells2(i).frame(end);
        else
            cdata(i).dsize = [];
            cdata(i).td = [];
        end
        for j=1:length(cells(i).frame)
            obj = frame3(cells(i).frame(j)).object(cells(i).object(j));
            cdata(i).area = [cdata(i).area obj.area*pixel2micron^2];
            cdata(i).clen = [cdata(i).clen obj.perim*pixel2micron];
            % SA & V
            cdata(i).SA = [cdata(i).SA ...
                surface_area_from_mesh(obj.pill_mesh,obj.centerline)*pixel2micron^2];
            cdata(i).V = [cdata(i).V ...
                volume_from_mesh(obj.pill_mesh,obj.centerline)*pixel2micron^3];
        end
    end
end