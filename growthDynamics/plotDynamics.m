%% growth kinetics
close all
t = (1:fN).*dT;

% growth rate
gr_map = [];
gr_data = cell(1,fN);
bins = 0:0.003:0.042;
for ii=1:length(fData)
    fdata = fData{ii};
    for i=1:length(fdata)
        gr_map(:,i,ii) = histcounts(fdata(i).gr,bins);
        gr_data{i} = [gr_data{i} fdata(i).gr];
    end
end
gr_map_tot = sum(gr_map,3);
for i=1:fN
    I = (gr_data{i}>bins(1)) & (gr_data{i}<bins(end));
    gr_mean(i) = mean(gr_data{i}(I));
    gr_std(i) = std(gr_data{i}(I));
end
clear I

n = 3; % smoothin
gr_mean = smooth(gr_mean,n);
gr_std = smooth(gr_std,n);

figure(1), hold on
h=imagesc(max(gr_map_tot(:))-gr_map_tot);
box on
set(gca,'Ydir','normal')
ylim([0,max(bins)])
xlim([0,t(end)-n0-1])
set(h,'YData',bins);
set(h,'XData',t);
colormap gray
plot(t,gr_mean,'linewidth',2,'color',[0.64,0.08,0.18]);
plot(t,gr_mean+gr_std,'lineWidth',0.7,'color',[0.77,0.1,0.22],'lineStyle','--')
plot(t,gr_mean-gr_std,'lineWidth',0.7,'color',[0.77,0.1,0.22],'lineStyle','--')

xlabel('time (min)')
ylabel('growth rate (dA/dt)/A')

clear n

%% growth rate vs. cell area over frames (growth kinetics)
close all

all_data = [];
for ii=1:length(fData)
    fdata = fData{ii};
    fnum = [];
    for i=1:length(fdata)
        fnum = [fnum i.*ones(1,length(fdata(i).gr))];
    end
    % collect all data points
    all_data = [all_data; cat(2,fdata.A)', cat(2,fdata.gr)', cat(2,fdata.dAdt)', fnum'];
end
all_data(:,4) = all_data(:,4).*dT;

figure(2)
scatter(all_data(:,1),all_data(:,2),0.7,all_data(:,4),'MarkerFaceAlpha',0.5)
colormap('copper')
colorbar
ylim([-0.01,0.05])
xlabel('cell area ($\mu m^2$)','interpreter','latex')
ylabel('(dA/dt)/A')
box on
figure(3)
scatter(all_data(:,1),all_data(:,3),0.7,all_data(:,4),'MarkerFaceAlpha',0.5)
colormap('copper')
colorbar
ylim([-0.05,0.16])
xlabel('cell area ($\mu m^2$)','interpreter','latex')
ylabel('dA/dt')
box on
% binning (every ? frames)
A_bin = linspace(min(all_data(:,1)),max(all_data(:,1)),9);
T_bin = 0:50:150;
binning_data = zeros(length(T_bin)-1,length(A_bin)-1);
binning_sem = size(binning_data);
for i=1:length(T_bin)-1
    temp = all_data( all_data(:,4)>T_bin(i) & all_data(:,4)<=T_bin(i+1),:);
    temp((temp(:,2)<-0.02)|(temp(:,2)>0.06),:) = [];
    for j=1:length(A_bin)-1
        I = temp(:,1)>A_bin(j) & temp(:,1)<=A_bin(j+1);
        binning_data(i,j) = mean(temp(I,2));
        binning_sem(i,j) = std(temp(I,2))/sqrt(length(temp(I,2)));
    end
end
figure(2), hold on
for i=1:length(T_bin)-1
    %plot(A_bin(1:end-1)+0.5*(A_bin(2)-A_bin(1)),binning_data(i,:))
    errorbar(A_bin(1:end-1)+0.5*(A_bin(2)-A_bin(1)),binning_data(i,:),binning_sem(i,:),...
        'linewidth',2)
end
clear all_data A_bin T_bin binning_data binning_sem temp I

%% single cell growth trajectories (growth kinetics)
close all

figure(4), hold on
cmap = colormap(copper(fN));
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if length(cdata(j).frame)==length(cdata(j).area)
            if ~isempty(cdata(j).tb)
                plot((cdata(j).frame-cdata(j).frame(1)).*dT,cdata(j).area,'color',cmap(cdata(j).tb,:))
            else
                plot((cdata(j).frame-cdata(j).frame(1)).*dT,cdata(j).area,'color',cmap(1,:))
            end
        end
    end
end
xlabel('time since birth (min)')
ylabel('cell area ($\mu m^2$)','interpreter','latex')
colormap('copper')
caxis([1,fN])
%caxis([1,30])
colorbar

figure(5), hold on % log scale
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if length(cdata(j).frame)==length(cdata(j).area)
            if ~isempty(cdata(j).tb)
                plot((cdata(j).frame-cdata(j).frame(1)).*dT,cdata(j).area,'color',cmap(cdata(j).tb,:))
            else
                plot((cdata(j).frame-cdata(j).frame(1)).*dT,cdata(j).area,'color',cmap(1,:))
            end
        end
    end
end
set(gca, 'YScale', 'log')
xlabel('time since birth (min)')
ylabel('cell area ($\mu m^2$)','interpreter','latex')
colormap('copper')
caxis([1,fN])
%caxis([1,30])
colorbar

%% size control models
close all

b = []; d = []; ti = []; tb = [];
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if ~isempty(cdata(j).bsize) && ~isempty(cdata(j).dsize)
            b = [b cdata(j).bsize];
            d = [d cdata(j).dsize];
            ti = [ti cdata(j).td - cdata(j).tb];
            tb = [tb cdata(j).tb];
        end
    end
end
tb = tb.*dT;
ti = ti.*dT;

% birth size vs. division size
figure(6), hold on
scatter(b,d,10,tb,'filled');
xlabel('birth size ($\mu m^2$)','interpreter','latex');
ylabel('division size ($\mu m^2$)','interpreter','latex');
colormap('copper');
colorbar
axis equal, box on
%caxis([1,30])
% binning (every ? frames)
b_bin = [0.8,1.2,1.6,2.0,2.4,2.8];
T_bin = [0,150];
binning_data = zeros(length(T_bin)-1,length(b_bin)-1);
binning_sem = size(binning_data);
for i=1:length(T_bin)-1
    temp = [b( tb>=T_bin(i) & tb<T_bin(i+1))' d( tb>=T_bin(i) & tb<T_bin(i+1))'];
    for j=1:length(b_bin)-1
        I = temp(:,1)>b_bin(j) & temp(:,1)<=b_bin(j+1);
        binning_data(i,j) = mean(temp(I,2));
        binning_sem(i,j) = std(temp(I,2))/sqrt(length(temp(I,2)));
    end
end
for i=1:length(T_bin)-1
    errorbar(b_bin(1:end-1)+0.5*(b_bin(2)-b_bin(1)),binning_data(i,:),binning_sem(i,:),...
        'linewidth',2)
end
xlim([0.7,3.2])
%clear b_bin T_bin binning_data binning_sem temp I

% birth size vs. size increase
figure(61)
scatter(b,d-b,30,tb,'filled');
xlabel('birth size ($\mu m^2$)','interpreter','latex');
ylabel('increased size ($\mu m^2$)','interpreter','latex');
colormap('copper');
colorbar
axis equal, box on
%caxis([1,30])

figure(7), hold on % later time points
I = tb>100;
scatter(b(I),d(I),30,tb(I),'filled');
% p = polyfit(b(I),d(I),1);
% x = [min(b(I)),max(b(I))];
% plot(x,p(1).*x + p(2))
xlabel('birth size ($\mu m^2$)','interpreter','latex')
ylabel('division size ($\mu m^2$)','interpreter','latex')
colormap('copper')
colorbar
caxis([min(tb),max(tb)])
axis equal, box on

figure(71), hold on % earlier time points
I = tb<50;
scatter(b(I),d(I),30,tb(I),'filled');
% p = polyfit(b(I),d(I),1);
% x = [min(b(I)),max(b(I))];
% plot(x,p(1).*x + p(2))
xlabel('birth size ($\mu m^2$)','interpreter','latex')
ylabel('division size ($\mu m^2$)','interpreter','latex')
colormap('copper')
colorbar
caxis([min(tb(I)),max(tb(I))])
axis equal, box on

% birth size vs. inter-division time
figure(8)
scatter(b,ti,30,tb,'filled')
colormap('copper')
colorbar
xlabel('birth size ($\mu m^2$)','interpreter','latex')
ylabel('inter-division time (min)')
box on
%caxis([1,30])

%% shape control
close all

% single cell contour length
figure(9), hold on
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        plot((1:length(cdata(j).clen)).*dT,cdata(j).clen)
    end
end
xlabel('time since birth (min)')
ylabel('contour length $(\mu m)$','interpreter','latex')

% single cell SA
figure(10), hold on
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        plot((1:length(cdata(j).SA)).*dT,cdata(j).SA)
    end
end
xlabel('time since birth (min)')
ylabel('surface area $(\mu m^2)$','interpreter','latex')

% single cell V
figure(11), hold on
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if ~isempty(cdata(j).tb)
            plot(cdata(j).V,'color',cmap(cdata(j).tb,:))
        else
            plot(cdata(j).V,'color',cmap(1,:))
        end
    end
end
xlabel('time since birth (min)')
ylabel('cell volume $(um^3)$','interpreter','latex')
colormap('copper')
caxis([1,fN])
colorbar

figure(12), hold on % log scale
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if ~isempty(cdata(j).tb)
            plot(cdata(j).V,'color',cmap(cdata(j).tb,:))
        else
            plot(cdata(j).V,'color',cmap(1,:))
        end
    end
end
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('cell volume (um^3)')
colormap('copper')
caxis([1,fN])
colorbar

%% average shape parameters (area, contour length, SA, V)
close all

mean_area = [];
mean_clen = [];
mean_SA = [];
mean_V = [];
std_area = [];
std_clen = [];
for i=1:fN
    temp1 = []; temp2 = []; temp3 = []; temp4 = []; temp5 = [];
    for ii=1:length(fData)
        if i<=length(fData{ii})
            temp1 = [temp1 fData{ii}(i).area];
            temp2 = [temp2 fData{ii}(i).clen];
            temp3 = [temp3 fData{ii}(i).SA];
            temp4 = [temp4 fData{ii}(i).V];
            temp5 = [temp5 fData{ii}(i).SA./fData{ii}(i).V];
        end
    end
    mean_area(i) = mean(temp1);
    mean_clen(i) = mean(temp2);
    mean_SA(i) = mean(temp3);
    mean_V(i) = mean(temp4);
    std_area(i) = std(temp1);
    std_clen(i) = std(temp2);
end
figure(14), hold on
fill([(1:length(mean_area)).*dT, fliplr((1:length(mean_area)).*dT)],...
    [mean_area+std_area fliplr(mean_area-std_area)],[0.9,0.9,0.9],...
    'LineStyle','none')
plot((1:length(mean_area)).*dT,mean_area,'k');
xlabel('time (min)')
ylabel('mean cell area ($\mu m^2$)','interpreter','latex')
box on

%% SA/V dynamics
close all

figure(15), hold on
plot(mean_SA)
xlabel('frame')
ylabel('mean surface area (um^2)')

figure(16), hold on
plot(mean_V)
xlabel('frame')
ylabel('mean volume (um^3)')

% single cell SA/V
figure(17), hold on
cmap = colormap(copper(fN));
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if length(cdata(j).frame)==length(cdata(j).SA) && length(cdata(j).SA)>6
            %plot(cdata(j).frame,cdata(j).SA./cdata(j).V,'g')
            %plot(cdata(j).frame,cdata(j).SA./cdata(j).V,'color',cmap(cdata(j).frame(1),:))
            plot(cdata(j).frame,cdata(j).SA./cdata(j).V)
%             if ~isempty(cdata(j).tb)
%                 plot(cdata(j).SA./cdata(j).V,'color',cmap(cdata(j).tb,:))
%             else
%                 plot(cdata(j).SA./cdata(j).V,'color',cmap(1,:))
%             end
        end
    end
end
load('colorMap.mat')
plot(mean_SA./mean_V,'color',colorMap(1,:),'LineWidth',2)
xlabel('frame')
ylabel('SA/V')
% colorbar
% caxis([1,fN])

%% initial cells
% cell area vs. frame
load('colorMap.mat')
figure(50), hold on % log scale
for ii=1:length(cData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if cdata(j).frame(1)<15 && length(cdata(j).frame)==length(cdata(j).area)
            plot(cdata(j).frame,cdata(j).area,'color',colorMap(ii,:))
        end
    end
end
set(gca, 'YScale', 'log')
xlabel('time since birth')
ylabel('cell area (um^2)')

figure(51), hold on
for ii=1:length(fData)
    cdata = cData{ii};
    for j=1:length(cdata)
        if cdata(j).frame(1)<=cdata(end).frame(1)+10 && isempty(cdata(j).tb)...
                && length(cdata(j).frame)==length(cdata(j).area)
            plot(cdata(j).frame,cdata(j).area)
        end
    end
end
set(gca, 'YScale', 'log')
xlabel('frame')
ylabel('cell area (um^2)')
title('initial cells')