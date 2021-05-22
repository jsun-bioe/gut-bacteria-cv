function pointed = traceBack(framePointer,objectPointer,frame,f)
% Trace back to a previous frame (usually 2 frames back) to find out cells 
% that have already been pointed at in the tree of the current cell
% Note that this is not a perfect algorithm, we can mistakenly reject 
% certain links in extreme cases, e.g. A->(B,missing)->C (Can it be fixed?)!

f0 = framePointer;
obj0 = objectPointer;
flag = 0;

while flag==0
    % if the pointers are empty (no mother cell yet), return []
    if isempty(f0)
        pointed = [];
        break
    end
    
    f1 = [];
    obj1 = [];
    for j=1:length(f0)
        if f0(j)>f
            f1 = [f1 frame(f0(j)).object(obj0(j)).framePointer];
            obj1 = [obj1 frame(f0(j)).object(obj0(j)).objectPointer];
        elseif f0(j)==f
            f1 = [f1 f0(j)];
            obj1 = [obj1 obj0(j)];
        end
    end
    I = f1<f;
    f1(I) = [];
    obj1(I) = [];
    f0 = f1;
    obj0 = obj1;
    if max(f0)==f
        flag = 1;
    end
end
pointed = unique(obj0);

end

