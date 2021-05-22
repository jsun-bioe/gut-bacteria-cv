function visualizeTree(cells,ISFIXED_)
% Visualize tree structure (raw or fixed) obtained in previous steps

load('colorMap.mat');

p = [];
for i=1:length(cells)
    mother = cells(i).mother;
    if isempty(mother)
        p(i) = 0;
    elseif length(mother)>1 %merging
        p(i) = 0; %treat as a starter (will be corrected in later steps)
        %p(i) = mother(1);
    else
        p(i) = mother;
    end
end

% plot tree structure

% figure, hold on
% treeplot(p)
[x,y] = treelayout(p);
% for i=1:length(x)
%     text(x(i),y(i)-0.05,num2str(length(cells)-i+1))
% end

figure, hold on
for i=1:length(cells)
    if ISFIXED_ && abs(cells(i).connected)==1
        continue
    end
    mother = cells(i).mother;
%     if length(mother)>1 % merging
%         x(i) = mean(x(mother));
%     end
    scatter(cells(i).frame,x(i).*ones(1,length(cells(i).frame)),5,colorMap(1,:),'filled')
    plot([cells(i).frame(1),cells(i).frame(end)],[x(i),x(i)])
    for j=1:length(mother)
        plot([cells(mother(j)).frame(end),cells(i).frame(1)],[x(mother(j)),x(i)])
    end
end
set(gca,'ytick',[])
set(gca,'ycolor','none')
xlabel('frame')

end

